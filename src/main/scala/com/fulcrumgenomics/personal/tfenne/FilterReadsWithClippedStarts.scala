/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.personal.tfenne

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.{SAMFileWriterFactory, SAMRecord, SamReaderFactory}

@clp(group=ClpGroups.Personal, description=
  """
    |Filters out all reads from templates where either read 1 or read 2 has >= N bases clipped from the start.
    |
    |If the input is queryname sorted or grouped the output file will be in the same order as the input file.
    |Otherwise the output file will be queryname sorted.
  """)
class FilterReadsWithClippedStarts
( @arg(flag="i", doc="Input aligned SAM or BAM file.") val input: PathToBam,
  @arg(flag="o", doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag="m", doc="Maximum clipping at the start of reads.") val maxClipping: Int = 5,
  @arg(flag="t", doc="Temporary directory to use when sorting.") val tmpDir: Option[DirPath] = None

) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in      = SamReaderFactory.make.open(input)
    val sorting = in.getFileHeader.getSortOrder != SortOrder.queryname && in.getFileHeader.getGroupOrder != GroupOrder.query
    val out     = new SAMFileWriterFactory().setCreateIndex(true).makeWriter(in.getFileHeader, !sorting, output.toFile, null)

    var (total, retained) = (0L, 0L)
    val progress          = new ProgressLogger(logger, unit=5e6.toInt)

    Bams.templateIterator(in, tmpDir=tmpDir).foreach { template =>
      // Only look at templates that actually have a primary read!
      if (template.r1.isDefined || template.r2.isDefined) {
        total += template.readCount

        if (!rejectRead(template.r1, maxClipping) && !rejectRead(template.r2, maxClipping)) {
          retained += template.readCount
          template.allReads.foreach(out.addAlignment)
        }
      }

      template.allReads.foreach(progress.record)
    }

    out.close()
    in.safelyClose()
    val fraction = retained / total.toDouble
    logger.info(f"Retained ${retained}%,d of ${total}%,d total reads (${fraction * 100}%2.3f)")
  }

  /** Reject a read if the amount of clipping at the start of the read is greater than maxClipping. */
  def rejectRead(rec: Option[SAMRecord], maxClipping: Int): Boolean = {
    rec match {
      case None    => false
      case Some(r) =>
        val cigar = {
          val elems = r.getCigar.getCigarElements.toSeq
          if (r.getReadNegativeStrandFlag) elems.reverse else elems
        }

        val startClipping = cigar.takeWhile(_.getOperator.isClipping).map(_.getLength).sum
        startClipping > maxClipping
    }
  }
}
