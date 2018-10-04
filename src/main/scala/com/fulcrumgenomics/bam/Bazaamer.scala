/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef.{PathToBam, PathToFastq, SafelyClosable, javaIterableToIterator, javaIteratorAsScalaIterator}
import com.fulcrumgenomics.bam.api._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fastq.{FastqRecord, FastqWriter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.samtools.util.{CoordMath, Interval, IntervalList, SequenceUtil}
import htsjdk.samtools.{BAMFileReader, QueryInterval, SAMRecord}

import scala.collection.mutable

object Region {
  def apply(region: String): Region = {
    region.split(':') match {
      case Array(contig, rest) =>
        val Array(start, end) = rest.split('-').map(_.toInt)
        Region(contig, start, end)
      case _ => throw new IllegalArgumentException(s"Could not parser region: '$region'.")
    }
  }
}
case class Region(contig: String, start: Int, end: Int) {
  def overlaps(contig: String, start: Int, end: Int): Boolean = {
    if (contig != this.contig) false
    else CoordMath.overlaps(start, end, this.start, this.end)
  }
  def overlaps(other: Region): Boolean = this.overlaps(other.contig, other.start, other.end)
}

@clp(group=ClpGroups.SamOrBam, description =
  """
    |Bazaam
    |
    |Extracts all reads that overlap the given region.  Also extracts the mates if they map outside the given region.
  """)
class Bazaamer
( @arg(flag='i', doc="Input BAM.") val input: PathToBam,
  @arg(flag='o', doc="Output FASTQ") val output: PathToFastq,
  @arg(flag='r', doc="Region to extract `<chr>:<start>-<end>`") region: Region
) extends FgBioTool with LazyLogging {

  implicit class MateEnd(rec: SamRecord) {
    def toMateEnd: Int = rec.mateEnd.getOrElse { throw new IllegalStateException("Cannot get mate end: MC tag required.")}
  }

  private def toFastqRecord(rec: SamRecord): FastqRecord = {
    toFastqRecord(rec.name, rec.positiveStrand, rec.basesString, rec.qualsString, rec.paired, rec.firstOfPair)
  }

  private def toFastqRecord(name: String, positiveStrand: Boolean, basesString: String, qualsString: String,
                            paired: Boolean, firstOfPair: Boolean): FastqRecord = {
    FastqRecord(
      name = name,
      bases = if (positiveStrand) basesString else SequenceUtil.reverseComplement(basesString),
      quals = if (positiveStrand) qualsString else qualsString.reverse,
      readNumber = if (!paired) None else if (firstOfPair) Some(1) else Some(2)
    )
  }

  override def execute(): Unit = {
    val reader = SamSource(input)
    require(reader.indexed, "Input BAM was not indexed")

    val mateIntervals = new mutable.HashSet[Interval]()
    val reads = new mutable.HashMap[String,FastqRecord]()
    val writer = FastqWriter(output)
    var readsWritten = 0L

    // A little method to write to the FASTQ if both reads are found, otherwise cache in `reads`
    def processSamRecord(rec: SamRecord): Unit = {
      reads.remove(rec.name) match {
        case None     => reads.put(rec.name, toFastqRecord(rec))
        case Some(fq) =>
          val (fq1, fq2) = if (rec.firstOfPair) (toFastqRecord(rec), fq) else (fq, toFastqRecord(rec))
          writer += fq1
          writer += fq2
          readsWritten += 2
      }
    }

    // Step 1: get all the reads in the specified region.  If the read is paired and its mate falls outside the region,
    // then we'll need to cache the read and find its mate later.
    val queryProgress = ProgressLogger(logger, verb = "query region")
    reader
      .query(region.contig, region.start, region.end, QueryType.Overlapping)
      .filterNot(r => r.secondary || r.supplementary)
      .foreach { rec =>
        // Check to see if the mate does not overlap the region being queried. If it doesn't overlap, we'll need to
        // retrieve it
        if (rec.paired && rec.mateStart != SamRecord.UnmappedStart && !region.overlaps(rec.mateRefName, rec.mateStart, rec.toMateEnd)) {
          mateIntervals += new Interval(rec.mateRefName, rec.mateStart, rec.toMateEnd)
        }
        processSamRecord(rec)
        queryProgress.record(rec)
      }
    logger.info(f"Wrote $readsWritten%,d reads with mates inside the region.")
    logger.info(f"${reads.size}%,d reads with mates outside of region.")

    // Step 2: Find the mates that fall outside the region.  Note: we merge all the intervals across all mates that fall
    // outside the originally specified region, to reduce the overall # of index queries.
    if (reads.nonEmpty) {
      logger.info(f"Found ${mateIntervals.size}%,d mate intervals")
      val queryIntervals = {
        val intervalList = new IntervalList(reader.dict)
        mateIntervals.foreach { i => intervalList.add(i) }
        intervalList.sorted.uniqued(false).getIntervals
      }
      logger.info(f"Found ${queryIntervals.length}%,d query intervals")
      queryIntervals.foreach { interval =>
        logger.info(f"    ${interval.toString}")
      }

      val overlappingProgress = ProgressLogger(logger, verb = "query mates", unit = 1000)
      reader.query(queryIntervals, QueryType.Contained).foreach { rec =>
        overlappingProgress.record(rec)
        if (reads.contains(rec.name)) processSamRecord(rec)
      }
    }
    logger.info(f"Wrote $readsWritten%,d reads total.")
    require(reads.isEmpty, f"Did not find mates for ${reads.size}%,d reads")

    reader.safelyClose()
    writer.close()
  }
}
