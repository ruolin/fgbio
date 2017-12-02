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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.commons.util.{CaptureSystemStreams, LogLevel, Logger}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import org.scalatest.OptionValues

class FlagStatTest extends UnitSpec with OptionValues with CaptureSystemStreams {


  private val builder: SamBuilder = new SamBuilder()
  private def metric: FlagStatMetric = FlagStatMetric(pf=true)

  "FlagStat.updateMetric" should "count stats for unmapped fragment reads" in {
    val m = metric
    val r = builder.addFrag(unmapped=true).value
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1)
  }

  it should "count stats for mapped fragment reads" in {
    val m = metric
    val r = builder.addFrag(start=1).value
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1)
  }

  it should "count stats for mapped secondary fragment reads" in {
    val m = metric
    val r = builder.addFrag(start=1).map(r => {r.secondary = true; r}).value
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, secondary=1)
  }

  it should "count stats for mapped supplementary fragment reads" in {
    val m = metric
    val r = builder.addFrag(start=1).map(r => {r.supplementary = true; r}).value
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, supplementary=1)
  }

  it should "count stats for mapped duplicate fragment reads" in {
    val m = metric
    val r = builder.addFrag(start=1).map(r => {r.duplicate = true; r}).value
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, duplicate=1)
  }

  it should "count stats for unmapped paired reads" in {
    val recs = builder.addPair(start1=1, start2=1, unmapped1=true, unmapped2=true)

    // first of pair
    {
      val m = metric
      FlagStat.updateMetric(m, recs.head)
      m shouldBe FlagStatMetric(pf=true, reads=1, paired=1, read_one=1)
    }

    // second of pair
    {
      val m = metric
      FlagStat.updateMetric(m, recs.last)
      m shouldBe FlagStatMetric(pf=true, reads=1, paired=1, read_two=1)
    }
  }

  it should "count stats for mapped paired reads" in {
    val m = metric
    val r = builder.addPair(start1=1, start2=1).head
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, paired=1, mapped_in_pair=1, read_one=1)
  }

  it should "count stats for mapped secondary paired reads" in {
    val m = metric
    val r = builder.addPair(start1=1, start2=1).head
    r.secondary = true
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, secondary=1)
  }

  it should "count stats for mapped supplementary paired reads" in {
    val m = metric
    val r = builder.addPair(start1=1, start2=1).head
    r.supplementary = true
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, supplementary=1)
  }

  it should "count stats for mapped duplicate paired reads" in {
    val m = metric
    val r = builder.addPair(start1=1, start2=1).head
    r.duplicate = true
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, paired=1, mapped_in_pair=1, read_one=1, duplicate=1)
  }

  it should "count stats for mapped proper paired reads" in {
    val m = metric
    val r = builder.addPair(start1=1, start2=1).head
    r.properlyPaired = true
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, paired=1, mapped_in_pair=1, proper_paired=1, read_one=1)
  }

  it should "count stats for singletons (mapped and paired read with mate unmapped)" in {
    val m = metric
    val r = builder.addPair(start1=1, start2=1, unmapped2=true).head
    r.mateUnmapped = true
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, paired=1, mapped_in_pair=1, read_one=1, singleton=1)
  }

  it should "count stats for inter-chromosomal mapped read pairs" in {
    val m = metric
    val r = builder.addPair(start1=1, start2=1).head
    r.refIndex = 1
    r.mapq = 4
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, paired=1, mapped_in_pair=1, read_one=1, inter_chrom=1)
  }

  it should "count stats for inter-chromosomal mapped read pairs with mapping greater than five" in {
    val m = metric
    val r = builder.addPair(start1=1, start2=1).head
    r.refIndex = 1
    r.mapq = 5
    FlagStat.updateMetric(m, r)
    m shouldBe FlagStatMetric(pf=true, reads=1, mapped=1, paired=1, mapped_in_pair=1, read_one=1, inter_chrom=1, inter_chrom_mapq5=1)
  }

  "FlagStat.formatInSamtools" should "print output in SAMTools flagstat format" in {
    // NB: uses the records created in the test above!

    // all records are not pf
    {
      val metrics = Seq(true, false).map(pf => pf -> FlagStatMetric(pf=pf)).toMap
      builder.foreach { r =>
        r.pf = false
        FlagStat.updateMetric(metrics(r.pf), r)
        r.pf = true // reset to true
      }
      val actual = FlagStat.formatInSamtools(metrics)
      val expected =
        """0 + 23 in total (QC-passed reads + QC-failed reads)
          |0 + 2 secondary
          |0 + 2 supplementary
          |0 + 2 duplicates
          |0 + 23 mapped (N/A : 82.61%)
          |0 + 16 paired in sequencing
          |0 + 7 read1
          |0 + 9 read2
          |0 + 1 properly paired (N/A : 6.25%)
          |0 + 13 with itself and mate mapped
          |0 + 1 singletons (N/A : 6.25%)
          |0 + 2 with mate mapped to a different chr
          |0 + 1 with mate mapped to a different chr (mapQ>=5)
          |""".stripMargin
      actual shouldBe expected
    }

    // all records are pf
    {
      val metrics = Seq(true, false).map(pf => pf -> FlagStatMetric(pf=pf)).toMap
      builder.foreach { r =>
        r.pf = true
        FlagStat.updateMetric(metrics(r.pf), r)
      }
      val actual = FlagStat.formatInSamtools(metrics)
      val expected =
        """23 + 0 in total (QC-passed reads + QC-failed reads)
          |2 + 0 secondary
          |2 + 0 supplementary
          |2 + 0 duplicates
          |23 + 0 mapped (82.61% : N/A)
          |16 + 0 paired in sequencing
          |7 + 0 read1
          |9 + 0 read2
          |1 + 0 properly paired (6.25% : N/A)
          |13 + 0 with itself and mate mapped
          |1 + 0 singletons (6.25% : N/A)
          |2 + 0 with mate mapped to a different chr
          |1 + 0 with mate mapped to a different chr (mapQ>=5)
          |""".stripMargin
      actual shouldBe expected
    }

    // half are pf and half are not pf
    {
      val metrics = Seq(true, false).map(pf => pf -> FlagStatMetric(pf=pf)).toMap
      builder.foreach { r =>
        r.pf = false
        FlagStat.updateMetric(metrics(r.pf), r)
        r.pf = true
        FlagStat.updateMetric(metrics(r.pf), r)
      }
      val actual = FlagStat.formatInSamtools(metrics)
      val expected =
        """23 + 23 in total (QC-passed reads + QC-failed reads)
          |2 + 2 secondary
          |2 + 2 supplementary
          |2 + 2 duplicates
          |23 + 23 mapped (82.61% : 82.61%)
          |16 + 16 paired in sequencing
          |7 + 7 read1
          |9 + 9 read2
          |1 + 1 properly paired (6.25% : 6.25%)
          |13 + 13 with itself and mate mapped
          |1 + 1 singletons (6.25% : 6.25%)
          |2 + 2 with mate mapped to a different chr
          |1 + 1 with mate mapped to a different chr (mapQ>=5)
          |""".stripMargin
      actual shouldBe expected
    }
  }

  "FlagStat" should "run end-to-end" in {
    val input = builder.toTempFile()
    val output = makeTempFile("flagstat.", ".txt")

    val previousLogLevel = Logger.level
    Logger.level = LogLevel.Warning
    val stdout = captureStdout { () =>
      new FlagStat(input=input, output=Some(output)).execute()
    }
    Logger.level = previousLogLevel

    val expected =
      """23 + 0 in total (QC-passed reads + QC-failed reads)
        |2 + 0 secondary
        |2 + 0 supplementary
        |2 + 0 duplicates
        |23 + 0 mapped (82.61% : N/A)
        |16 + 0 paired in sequencing
        |7 + 0 read1
        |9 + 0 read2
        |1 + 0 properly paired (6.25% : N/A)
        |13 + 0 with itself and mate mapped
        |1 + 0 singletons (6.25% : N/A)
        |2 + 0 with mate mapped to a different chr
        |1 + 0 with mate mapped to a different chr (mapQ>=5)
        |""".stripMargin
    stdout shouldBe expected
  }
}
