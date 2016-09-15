/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.metagenomics

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools._
import htsjdk.samtools.util.{SortingCollection, TrimmingUtil}

import scala.collection.JavaConversions.iterableAsScalaIterable


@clp(group=ClpGroups.Metagenomics, description=
"""
  |Takes one or more fastq files, and prepares a single fasta for kraken[1]. Individual reads
  |are quality trimmed using the phred/bwa method, and then any remaining low quality bases
  |are masked to Ns.  Read 1s and read 2s are the joined with a block of Ns between them
  |(of length kmerLength) so that the pair can be fed to kraken as a single read.
  |
  |Reads shorter than minReadLength after trimming are excluded from the output. For paired reads
  |if either read 1 or read 2 of a pair is shorter than minReadLength the pair is excluded.
  |
  |Input fastq files may be gzipped. Read 1 file(s) are required, read 2 file(s) are optional
  |but if provided must be equal in number to the read 1 files.
  |
  |[1] Kraken: http://ccb.jhu.edu/software/kraken/
""")
class PrepareKrakenInput
(
  @arg(flag="i", doc="Input BAM file.")                         val input: PathToBam,
  @arg(flag="o", doc="Output (singular) fasta file to write.")  val output: PathToFasta = Io.StdOut,
  @arg(flag="q", doc="Trim/mask bases below quality q.")        val minQuality: Int = 20,
  @arg(flag="k", doc="Kmer size of kraken database.")           val kmerSize: Int = 31,
  @arg(flag="m", doc="Discard template with fewer than this kmers after processing.") val minKmers: Int = 50
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  // String of N bases that is used between reads for the sample template in the output
  private val spacer = "N" * kmerSize

  override def execute(): Unit = {
    val out = Io.toWriter(output)

    // Utility method to write a fasta record to an output writer
    def writeFasta(name: String, bases: String): Unit = {
      out.write('>')
      out.write(name)
      out.newLine()
      out.write(bases)
      out.newLine()
    }

    val iterator = byQueryIterator(input).filterNot(_.isSecondaryOrSupplementary).bufferBetter
    val progress = new ProgressLogger(logger)
    while (iterator.hasNext) {
      val name = iterator.head.getReadName
      val baseString = iterator.takeWhile(_.getReadName == name).toSeq
        .map(r => trimAndMask(r.getReadBases, r.getBaseQualities))
        .map(bs => new String(bs))
        .mkString(spacer)

      if (hasMinKmers(baseString)) writeFasta(name, baseString)
      progress.record()
    }

    out.close()
  }

  /** Returns an iterator over records in queryname order, sorting if necessary. */
  private def byQueryIterator(bam: PathToBam): Iterator[SAMRecord] = {
    val in = SamReaderFactory.make().open(bam)
    val header = in.getFileHeader

    if (header.getSortOrder == SortOrder.queryname || header.getGroupOrder == GroupOrder.query) {
      in.toIterator
    }
    else {
      logger.info("Sorting records into queryname order.")
      val progress = new ProgressLogger(logger, verb="Sorted")
      val sorter = SortingCollection.newInstance[SAMRecord](classOf[SAMRecord], new BAMRecordCodec(header), new SAMRecordQueryNameComparator, 1e6.toInt)
      in.foreach { rec =>
        shrinkRead(rec)
        sorter.add(rec)
        progress.record(rec)
      }
      logger.info("Done sorting records into queryname order.")
      sorter.toIterator
    }
  }

  /** Does what it can to skinny the read down as much as possible before sorting. */
  private def shrinkRead(rec: SAMRecord): Unit = {
    SAMUtils.makeReadUnmapped(rec)
    val rg = rec.getAttribute("RG")
    rec.clearAttributes()
    rec.setAttribute("RG", rg)
  }

  /**
    * Trims the end of the read based on quality and then masks any remaining low quality bases to Ns.
    *
    * @param bases the read bases as an array of bytes
    * @param quals the base qualities as an array of bytes (integer phred scores)
    * @return a trimmed and masked version of the bases
    * */
  private def trimAndMask(bases: Array[Byte], quals: Array[Byte]): Array[Byte] = {
    val trimPoint = TrimmingUtil.findQualityTrimPoint(quals, this.minQuality)
    forloop (from=0, until=trimPoint) { i =>
        if (quals(i) < this.minQuality) bases(i) = 'N'
    }

    bases.take(trimPoint)
  }

  /** Counts how many kmers the read has that do not contain any Ns. */
  private def hasMinKmers(s: String): Boolean = {
    if (s.length < this.kmerSize) {
      false
    }
    else {
      var goodLength = 0
      var kmers = 0
      forloop (from=0, until=s.length) { i =>
        if (s.charAt(i) == 'N') goodLength = 0
        else goodLength += 1

        if (goodLength >= this.kmerSize) kmers += 1
      }

      kmers >= this.minKmers
    }
  }
}
