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

package com.fulcrumgenomics

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamSource, SamWriter}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.testing.UnitSpec
import enumeratum.EnumEntry
import htsjdk.samtools.{SAMFileHeader, SAMReadGroupRecord}

import scala.collection.immutable.IndexedSeq

class FgBioDefTest extends UnitSpec {
  private case class SampleAndLibrary(sample: String, library: String)

  private def makeEmptyBam(sampleAndLibraries: SampleAndLibrary*): PathToBam = {
    val path = makeTempFile("tmp", ".bam")
    val header     = {
      val h = new SAMFileHeader()
      sampleAndLibraries.zipWithIndex.foreach { case (sal, idx) =>
        val rg = new SAMReadGroupRecord(s"$idx")
        rg.setSample(sal.sample)
        rg.setLibrary(sal.library)
        h.addReadGroup(rg)
      }
      SamOrder.Unsorted.applyTo(h)
      h
    }

    val writer = SamWriter(path, header, sort=Some(SamOrder.Unsorted))
    writer.close()
    path
  }

  "plotDescription" should "return the sample and library name if they are unique in the read groups" in {
    // one read group
    {
      val path   = makeEmptyBam(SampleAndLibrary("S", "L"))
      val reader = SamSource(path)
      FgBioDef.plotDescription(reader, path) shouldBe "S / L"
      reader.close()
    }

    // three read groups
    {
      val path   = makeEmptyBam(SampleAndLibrary("S", "L"), SampleAndLibrary("S", "L"), SampleAndLibrary("S", "L"))
      val reader = SamSource(path)
      FgBioDef.plotDescription(reader, path) shouldBe "S / L"
      reader.close()
    }
  }

  it should "return the input name (extension removed) if no read groups are found" in {
    val path   = makeEmptyBam()
    val reader = SamSource(path)
    FgBioDef.plotDescription(reader, path) shouldBe PathUtil.basename(path, trimExt=true).toString
    reader.close()
  }

  it should "return the input name (extension removed) if sample name is not unique in the read groups" in {
    val path   = makeEmptyBam(SampleAndLibrary("S1", "L"), SampleAndLibrary("S2", "L"), SampleAndLibrary("S1", "L"))
    val reader = SamSource(path)
    FgBioDef.plotDescription(reader, path) shouldBe PathUtil.basename(path, trimExt=true).toString
    reader.close()
  }

  it should "return the input name (extension removed) if library name is not unique in the read groups" in {
    val path   = makeEmptyBam(SampleAndLibrary("S", "L1"), SampleAndLibrary("S", "L2"), SampleAndLibrary("S", "L1"))
    val reader = SamSource(path)
    FgBioDef.plotDescription(reader, path) shouldBe PathUtil.basename(path, trimExt=true).toString
    reader.close()
  }

  sealed trait CustomEnum extends EnumEntry
  object CustomEnum extends FgBioEnum[CustomEnum] {
    def values: IndexedSeq[CustomEnum] = findValues
    case object Foo extends CustomEnum
    case object Bar extends CustomEnum
  }

  "FgbioDef.FgBioEnum" should "return the correct enum by name" in {
    CustomEnum("Foo") shouldBe CustomEnum.Foo
    CustomEnum("foo") shouldBe CustomEnum.Foo
    CustomEnum("fOO") shouldBe CustomEnum.Foo
    CustomEnum("bar") shouldBe CustomEnum.Bar
  }

  "FgbioDef.GenomicRange()" should "parse a `<chr>:<start>-<end>` correctly" in {
    an[IllegalArgumentException] should be thrownBy FgBioDef.GenomicRange("chr1")
    an[IllegalArgumentException] should be thrownBy FgBioDef.GenomicRange("chr1-")
    an[IllegalArgumentException] should be thrownBy FgBioDef.GenomicRange("chr1-2")
    an[IllegalArgumentException] should be thrownBy FgBioDef.GenomicRange("chr1:")
    an[IllegalArgumentException] should be thrownBy FgBioDef.GenomicRange("chr1:1-")
    an[NumberFormatException] should be thrownBy FgBioDef.GenomicRange("chr1:1-A")
    an[NumberFormatException] should be thrownBy FgBioDef.GenomicRange("chr1:B-1")
    an[NumberFormatException] should be thrownBy FgBioDef.GenomicRange("chr1:2-,")

    FgBioDef.GenomicRange("chr1:1")   shouldBe GenomicRange("chr1", 1, 1)
    FgBioDef.GenomicRange("chr1:1-2") shouldBe GenomicRange("chr1", 1, 2)
    FgBioDef.GenomicRange("foo:11-5") shouldBe GenomicRange("foo", 11, 5)
    FgBioDef.GenomicRange("chr2:1,234-5,678") shouldBe GenomicRange("chr2", 1234, 5678)
    FgBioDef.GenomicRange("chr1::1") shouldBe GenomicRange("chr1:", 1, 1)
    FgBioDef.GenomicRange("foo:bar:11-5") shouldBe GenomicRange("foo:bar", 11, 5)
  }
}
