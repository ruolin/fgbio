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
 *
 */

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.util.ReadStructure
import dagr.sopt.cmdline.ValidationException
import htsjdk.samtools.{SAMFileWriterFactory, SAMRecord}

/**
  * Tests for ExtractUmisFromBam.
  */
class ExtractUmisFromBamTest extends UnitSpec {

  def newBam = makeTempFile("extract_umis_from_bam_test.", ".bam")

  def annotateRecordFragment: SAMRecord = {
    new SamRecordSetBuilder(readLength=100).addFrag(name="Frag", start=1).map {rec => rec.setReadString("A" * 100); rec }.get
  }

  "ExtractUmisFromBam.annotateRecord" should "should not annotate with no molecular barcodes" in {
    val frag = annotateRecordFragment
    ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("100T"), Seq.empty) shouldBe 'empty
    frag.getReadLength shouldBe 100
  }

  it should "throw an exception when fewer molecular barcodes in the read structure than molecular SAM tags are given then " in {
    val tags = Seq("AB", "CD")
    val frag = annotateRecordFragment
    an[Exception] should be thrownBy ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("10M80T"), tags)
  }

  it should "throw an exception when fewer molecular SAM tags are given than molecular barcodes in the read structure" in {
    val tags = Seq("AB", "CD")
    val frag = annotateRecordFragment
    an[Exception] should be thrownBy ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("10M10M10M70T"), tags)
  }

  it should "should annotate a single molecular barcode" in {
    val builder = new SamRecordSetBuilder(readLength=100)
    val frag = annotateRecordFragment
    ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("90T10M"), Seq("RX")) shouldBe "A" * 10
    frag.getReadLength shouldBe 90
  }

  it should "should annotate a two molecular barcode" in {
    val builder = new SamRecordSetBuilder(readLength=100)
    val frag = annotateRecordFragment
    ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("90T5M5M"), Seq("R1", "R2")) shouldBe "A" * 10
    frag.getReadLength shouldBe 90
  }

  it should "not accept a single read" in {
    val output = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    builder.addFrag(name="Frag1", start=1)
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }

  it should "not accept fragments" in {
    val output = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    builder.addFrag(name="Frag1", start=1)
    builder.addFrag(name="Frag2", start=1)
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }

  it should "not accept mapped reads" in {
    val output = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    builder.addPair(name="Pair", start1=1, start2=1)
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }

  it should "only not accept paired records with the first end before the second end" in {
    val input   = newBam
    val output  = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    builder.addPair(name="Pair", start1=1, start2=1)

    val writer = new SAMFileWriterFactory()
      .makeWriter(builder.header, true, input.toFile, null)
    builder.iterator.toSeq.reverse.foreach(writer.addAlignment)
    writer.close()

    an[Exception] should be thrownBy new ExtractUmisFromBam(input=input, output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }

  it should "not accept read pairs with the different names" in {
    val output  = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1)
    records.head.setReadName("PairNot")
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }

  it should "not accept multiple molecular barcode SAM tags and have a different number of molecular barcodes in the read structure" in {
    val output  = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    an[ValidationException] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T"), molecularBarcodeTags=Seq("A1", "B1", "C1")).execute()
  }

  it should "not accept duplicate molecular barcode SAM tags" in {
    val output  = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    an[ValidationException] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T"), molecularBarcodeTags=Seq("A1", "A1")).execute()
  }

  it should "extract the molecular barcodes and annotate them for a read pair in tags" in {
    val output  = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, record1Unmapped=true, record2Unmapped=true)
    records.foreach { rec => rec.setReadString("A" * 100) }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "1M99T"), molecularBarcodeTags=Seq("P1", "P2")).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      if (rec.getFirstOfPairFlag) {
        rec.getStringAttribute("P1") shouldBe "A" * 10
        rec.getReadLength shouldBe 90
        rec.getReadString shouldBe "A" * 90
      }
      else {
        rec.getStringAttribute("P2") shouldBe "A"
        rec.getReadLength shouldBe 99
        rec.getReadString shouldBe "A" * 99
      }
    }
  }

  it should "extract the molecular barcodes and annotate them for a read pair in tags and in the names" in {
    val output  = newBam
    val builder = new SamRecordSetBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, record1Unmapped=true, record2Unmapped=true)
    records.foreach { rec => rec.setReadString("A" * 100) }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "1M99T"), molecularBarcodeTags=Seq("P1", "P2"), annotateReadNames=true).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      if (rec.getFirstOfPairFlag) {
        rec.getStringAttribute("P1") shouldBe "A" * 10
        rec.getReadLength shouldBe 90
        rec.getReadString shouldBe "A" * 90
        rec.getReadName shouldBe "Pair+" + "A" * 11
      }
      else {
        rec.getStringAttribute("P2") shouldBe "A"
        rec.getReadLength shouldBe 99
        rec.getReadString shouldBe "A" * 99
        rec.getReadName shouldBe "Pair+" + "A" * 11
      }
    }
  }
}
