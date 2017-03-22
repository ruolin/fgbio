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
import com.fulcrumgenomics.util.{Io, Metric, NumericCounter, ProgressLogger}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SamReaderFactory

import scala.collection.mutable.ListBuffer

case class ClippingMetric(read: Int, end: String, clip_length: Int, count: Long) extends Metric

@clp(group=ClpGroups.Personal, description="Quantifies clipping at the starts and ends of reads.")
class CollectClippingMetrics
( @arg(flag="i", doc="Input BAM file.") val input: PathToBam,
  @arg(flag="m", doc="Minimum insert size.") val minInsertSize: Int = 160,
  @arg(flag="o", doc="Output metrics file.") val output: FilePath
) extends FgBioTool with LazyLogging {


  override def execute(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)

    val Seq(r1Start, r1End, r2Start, r2End) = Seq(1,2,3,4).map(_ => new NumericCounter[Int])
    val in = SamReaderFactory.make().open(input)
    val progress = new ProgressLogger(logger, unit=5e6.toInt)

    in.filter(Bams.isFrPair).filter(_.getInferredInsertSize >= minInsertSize).foreach {rec =>
      val cigar = {
        val elems = rec.getCigar.getCigarElements.toSeq
        if (rec.getReadNegativeStrandFlag) elems.reverse else elems
      }

      val (start, end) = if (rec.getFirstOfPairFlag) (r1Start, r1End) else (r2Start, r2End)
      val startClipping = cigar.takeWhile(_.getOperator.isClipping).map(_.getLength).sum
      val endClipping   = cigar.reverse.takeWhile(_.getOperator.isClipping).map(_.getLength).sum
      start.count(startClipping)
      end.count(endClipping)

      progress.record(rec)
    }

    val metrics = new ListBuffer[ClippingMetric]
    Seq((r1Start, 1, "start"), (r1End, 1, "end"), (r2Start, 2, "start"), (r2End, 2, "end")).foreach { case (counter, read, end) =>
        metrics ++= counter.map { case (clip, count) => ClippingMetric(read, end, clip, count) }
    }

    Metric.write(output, metrics)
  }
}
