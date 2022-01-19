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

package com.fulcrumgenomics.vcf.filtration

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec, VcfBuilder}
import com.fulcrumgenomics.vcf.api.VcfSource

class FilterSomaticVcfTest extends UnitSpec {
  // A pair of (paths to) VCFs for use in testing below
  private lazy val Seq(tumorOnlyVcf, tumorNormalVcf) = {
    Seq(Seq("tumor"), Seq("tumor", "normal")).map { samples =>
      val builder = VcfBuilder(samples)
      val includeNormal = samples.contains("normal")

      builder.add(pos=100, alleles=Seq("C", "A"),   gts=Seq(Gt("tumor", "C/A"),   Gt("normal", "C/C")).filter(g => g.sample == "tumor" || includeNormal))
      builder.add(pos=200, alleles=Seq("G", "A"),   gts=Seq(Gt("tumor", "G/A"),   Gt("normal", "G/G")).filter(g => g.sample == "tumor" || includeNormal))
      builder.add(pos=300, alleles=Seq("AAA", "A"), gts=Seq(Gt("tumor", "AAA/A"), Gt("normal", "AAA/AAA")).filter(g => g.sample == "tumor" || includeNormal))
      builder.add(pos=400, alleles=Seq("A", "T"),   gts=Seq(Gt("tumor", "A/T"),   Gt("normal", "A/A")).filter(g => g.sample == "tumor" || includeNormal))
      builder.add(pos=500, alleles=Seq("C", "G"),   gts=Seq(Gt("tumor", "C/G"),   Gt("normal", "C/C")).filter(g => g.sample == "tumor" || includeNormal))
      builder.toTempFile()
    }
  }

  // BAM file that has reads over each of the variant sites
  private lazy val bam = {
    val rlen    = 40
    val builder = new SamBuilder(readLength=rlen, baseQuality=40, sort=Some(SamOrder.Coordinate))

    /**
      * C>A SNP at 100 with signal for A-tailing and end repair artifacts
      * ATailingArtifactFilter should apply and output low p-value
      * EndRepairArtifactFilter should apply and output low p-value
      */
    for (pos <- 61 to 100; i <- 1 to 3) {
      builder.addPair(start1=pos,      start2=pos+rlen, bases1="C"*rlen, bases2="C"*rlen)
      builder.addPair(start1=pos-rlen, start2=pos     , bases1="C"*rlen, bases2="C"*rlen)
    }
    builder.addPair(start1=20, start2=62, bases1="A"*rlen, bases2="A"*rlen)
    builder.addPair(start1=21, start2=61, bases1="A"*rlen, bases2="A"*rlen)
    builder.addPair(start1=16, start2=63, bases1="A"*rlen, bases2="A"*rlen)

    /**
      * G>A SNP at 200 with low MAF but even reads
      * ATailingArtifactFilter should apply and output high p-value
      * EndRepairArtifactFilter should apply and output high p-value
      */
    for (pos <- 161 to 200; i <- 1 to 3) {
      builder.addPair(start1=pos,      start2=pos+rlen, bases1="G"*rlen, bases2="G"*rlen)
      builder.addPair(start1=pos-rlen, start2=pos     , bases1="G"*rlen, bases2="G"*rlen)
      if (i == 1 && pos % 4 == 0) {
        builder.addPair(start1=pos,      start2=pos+rlen, bases1="A"*rlen, bases2="A"*rlen)
        builder.addPair(start1=pos-rlen, start2=pos     , bases1="A"*rlen, bases2="A"*rlen)
      }
    }

    // No reads at the indel position because we won't look at those anyway

    /**
      * A>T SNP at 400 with signal for A-tailing and end repair artifacts
      * ATailingArtifactFilter should apply and output low p-value
      * EndRepairArtifactFilter should apply and output low p-value
      */
    for (pos <- 361 to 400; i <- 1 to 3) {
      builder.addPair(start1=pos,      start2=pos+rlen, bases1="A"*rlen, bases2="A"*rlen)
      builder.addPair(start1=pos-rlen, start2=pos     , bases1="A"*rlen, bases2="A"*rlen)
    }
    builder.addPair(start1=398, start2=430, bases1="T"*rlen, bases2="T"*rlen)
    builder.addPair(start1=398, start2=440, bases1="T"*rlen, bases2="T"*rlen)
    builder.addPair(start1=399, start2=437, bases1="T"*rlen, bases2="T"*rlen)
    builder.addPair(start1=400, start2=449, bases1="T"*rlen, bases2="T"*rlen)
    builder.addPair(start1=400, start2=441, bases1="T"*rlen, bases2="T"*rlen)

    /**
      * C>G SNP at 500
      * ATailingArtifactFilter should not apply
      * EndRepairArtifactFilter should apply and output low p-value
      */
    for (pos <- 461 to 500; i <- 1 to 3) {
      builder.addPair(start1=pos,      start2=pos+rlen, bases1="C"*rlen, bases2="C"*rlen)
      builder.addPair(start1=pos-rlen, start2=pos     , bases1="C"*rlen, bases2="C"*rlen)
    }
    builder.addPair(start1=498, start2=530, bases1="G"*rlen, bases2="G"*rlen)
    builder.addPair(start1=498, start2=540, bases1="G"*rlen, bases2="G"*rlen)
    builder.addPair(start1=499, start2=537, bases1="G"*rlen, bases2="G"*rlen)
    builder.addPair(start1=500, start2=549, bases1="G"*rlen, bases2="G"*rlen)
    builder.addPair(start1=500, start2=541, bases1="G"*rlen, bases2="G"*rlen)

    builder.toTempFile()
  }

  private val ATailInfoKey       = new ATailingArtifactLikelihoodFilter().readEndInfoLine.id
  private val ATailFilterKey     = new ATailingArtifactLikelihoodFilter().readEndFilterLine.id
  private val EndRepairInfoKey   = new EndRepairFillInArtifactLikelihoodFilter().readEndInfoLine.id
  private val EndRepairFilterKey = new EndRepairFillInArtifactLikelihoodFilter().readEndFilterLine.id

  "FilterSomaticVcf" should "work on an empty VCF" in {
    val emptyVcf    = VcfBuilder(Seq("tumor")).toTempFile()
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    new FilterSomaticVcf(input=emptyVcf, output=filteredVcf, bam=bam).execute()
    val reader = VcfSource(filteredVcf)
    reader.header.info.contains(ATailInfoKey) shouldBe true
    reader.header.filter.contains(ATailFilterKey) shouldBe true
    reader.header.info.contains(EndRepairInfoKey) shouldBe true
    reader.header.filter.contains(EndRepairFilterKey) shouldBe true
    reader.safelyClose()
  }

  it should "work on a single sample VCF" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    new FilterSomaticVcf(input=tumorOnlyVcf, output=filteredVcf, bam=bam).execute()
    val variants = readVcfRecs(filteredVcf)
    variants should have size 5
    variants(0).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(1).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(2).get[Float](ATailInfoKey).isDefined shouldBe false
    variants(3).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(4).get[Float](ATailInfoKey).isDefined shouldBe false
    variants.exists(_.filters.contains(ATailFilterKey)) shouldBe false // no threshold == no filtering

    variants(0).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(1).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(2).get[Float](EndRepairInfoKey).isDefined shouldBe false
    variants(3).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(4).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants.exists(_.filters.contains(EndRepairFilterKey)) shouldBe false // no threshold == no filtering
  }

  it should "fail on a single-sample VCF if an invalid sample name is provided" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    an[Exception] shouldBe thrownBy { new FilterSomaticVcf(input=tumorOnlyVcf, output=filteredVcf, bam=bam, sample=Some("WhoDis")).execute() }
  }

  it should "fail on a multi-sample VCF if no sample name is provided" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    an[Exception] shouldBe thrownBy { new FilterSomaticVcf(input=tumorNormalVcf, output=filteredVcf, bam=bam).execute() }
  }

  it should "fail on a multi-sample VCF if an invalid sample name is provided" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    an[Exception] shouldBe thrownBy { new FilterSomaticVcf(input=tumorNormalVcf, output=filteredVcf, bam=bam, sample=Some("WhoDis")).execute() }
  }

  it should "work on a multi-sample VCF if a sample name is given" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    new FilterSomaticVcf(input=tumorNormalVcf, output=filteredVcf, bam=bam, sample=Some("tumor")).execute()
    val variants = readVcfRecs(filteredVcf)
    variants should have size 5
    variants(0).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(1).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(2).get[Float](ATailInfoKey).isDefined shouldBe false
    variants(3).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(4).get[Float](ATailInfoKey).isDefined shouldBe false
    variants.exists(_.filters.contains(ATailFilterKey)) shouldBe false // no threshold == no filtering

    variants(0).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(1).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(2).get[Float](EndRepairInfoKey).isDefined shouldBe false
    variants(3).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(4).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants.exists(_.filters.contains(EndRepairFilterKey)) shouldBe false // no threshold == no filtering
  }

  it should "apply filters if filter-specific p-value thresholds are supplied" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    new FilterSomaticVcf(input=tumorOnlyVcf, output=filteredVcf, bam=bam, sample=Some("tumor"), aTailingDistance=Some(4), aTailingPValue=Some(0.001), endRepairFillInPValue = Some(0.001)).execute()
    val variants = readVcfRecs(filteredVcf)
    variants should have size 5
    variants(0).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(1).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(2).get[Float](ATailInfoKey).isDefined shouldBe false
    variants(3).get[Float](ATailInfoKey).isDefined shouldBe true
    variants(4).get[Float](ATailInfoKey).isDefined shouldBe false

    variants(0).filters.contains(ATailFilterKey) shouldBe true
    variants(1).filters.contains(ATailFilterKey) shouldBe false
    variants(2).filters.contains(ATailFilterKey) shouldBe false
    variants(3).filters.contains(ATailFilterKey) shouldBe true
    variants(4).filters.contains(ATailFilterKey) shouldBe false

    variants(0).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(1).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(2).get[Float](EndRepairInfoKey).isDefined shouldBe false
    variants(3).get[Float](EndRepairInfoKey).isDefined shouldBe true
    variants(4).get[Float](EndRepairInfoKey).isDefined shouldBe true

    variants(0).filters.contains(EndRepairFilterKey) shouldBe true
    variants(1).filters.contains(EndRepairFilterKey) shouldBe false
    variants(2).filters.contains(EndRepairFilterKey) shouldBe false
    variants(3).filters.contains(EndRepairFilterKey) shouldBe true
    variants(4).filters.contains(EndRepairFilterKey) shouldBe true
  }
}
