/**
  * Copyright (c) 2016, Fulcrum Genomics LLC
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *
  * 1. Redistributions of source code must retain the above copyright notice,
  * this list of conditions and the following disclaimer.
  *
  * 2. Redistributions in binary form must reproduce the above copyright notice,
  * this list of conditions and the following disclaimer in the documentation
  * and/or other materials provided with the distribution.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  * POSSIBILITY OF SUCH DAMAGE.
  */

package com.fulcrumgenomics.illumina

import com.fulcrumgenomics.util.Io
import org.scalatest.{FlatSpec, Matchers, OptionValues}

import java.nio.file.Paths

class SampleSheetTest extends FlatSpec with Matchers with OptionValues {
  private val testDir = Paths.get("src/test/resources/com/fulcrumgenomics/util/miseq/")

  "SampleSheet" should "load a simple sample sheet" in {
    val sampleSheet = SampleSheet(testDir.resolve("SampleSheet.csv"))
    sampleSheet.size shouldBe 12
    sampleSheet.foreach { sample =>
      sample.sampleId                    shouldBe "20000101-EXPID-" + sample.sampleOrdinal
      sample.sampleName                  shouldBe "Sample_Name_" + sample.sampleOrdinal
      sample.libraryId                   shouldBe "20000101-EXPID-" + sample.sampleOrdinal
      sample.project.value               shouldBe "Sample_Project_" + sample.sampleOrdinal
      sample.description.value           shouldBe "Description_" + sample.sampleOrdinal
      sample.i7IndexBases.value          shouldBe "GATTACAACGT"
      sample.i5IndexBases.value          shouldBe "GATTACA"
      sample.extendedAttributes.keySet should contain ("R1_Barcode_Bases".toUpperCase)
      sample.extendedAttributes.keySet should contain ("R2_Barcode_Bases".toUpperCase)
    }
  }

  it should "load a sample sheet with only required values" in {
    val sampleSheet = SampleSheet(testDir.resolve("SampleSheetOnlyRequired.csv"))
    sampleSheet.size shouldBe 12
    sampleSheet.foreach { sample =>
      sample.sampleId             shouldBe "20000101-EXPID-" + sample.sampleOrdinal
      sample.sampleName           shouldBe "Sample_Name_" + sample.sampleOrdinal
      sample.libraryId            shouldBe "20000101-EXPID-" + sample.sampleOrdinal
      sample.project.value        shouldBe "Sample_Project_" + sample.sampleOrdinal
      sample.description.value    shouldBe "Description_" + sample.sampleOrdinal
      sample.i7IndexBases         shouldBe None
      sample.i5IndexBases         shouldBe None
      sample.extendedAttributes.keySet should not contain "R1_Barcode_Bases"
      sample.extendedAttributes.keySet should not contain "R2_Barcode_Bases"
    }
  }

  it should "load a sample sheet and a subset of samples by lane" in {
    Seq(1, 2, 3).foreach { lane =>
      val sampleSheet = SampleSheet(testDir.resolve("SampleSheet.lanes.csv"), lane=Some(lane))
      sampleSheet.foreach { sample => sample.lane.value shouldBe lane }
      sampleSheet.size shouldBe 4
    }
    SampleSheet(testDir.resolve("SampleSheet.lanes.csv"), lane=Some(4)).isEmpty shouldBe true
    SampleSheet(testDir.resolve("SampleSheet.lanes.csv"), lane=None).size shouldBe 12
  }

  it should "load a sample sheet with empty lines at the end" in {
    val samples          = SampleSheet(testDir.resolve("SampleSheet.csv")).toSeq
    val sampleSheetLines = Io.readLines(testDir.resolve("SampleSheet.csv")).toSeq

    // Sample sheet with empty lines at the bottom
    val emptyLines = sampleSheetLines.iterator ++ Seq("", "")
    SampleSheet(lines = emptyLines, lane = None).toSeq should contain theSameElementsInOrderAs samples

    // Sample sheet with comma-only lines at the bottom (sometimes added by spreadsheet editors)
    val commaLines = sampleSheetLines.iterator ++ Seq(",,,,,,,", ",,,,,,,,")
    SampleSheet(lines = commaLines, lane = None).toSeq should contain theSameElementsInOrderAs samples

    // Sample sheet with mix of empty and comma-only lines at the bottom
    val commaAndEmptyLines = sampleSheetLines.iterator ++ Seq(",,,,,,,", ",,,,,,,,", "", "", ",,,,,,,,", "")
    SampleSheet(lines = commaAndEmptyLines, lane = None).toSeq should contain theSameElementsInOrderAs samples
  }

  it should "throw an exception when no sample section is found" in {
    an[Exception] should be thrownBy SampleSheet(testDir.resolve("SampleSheetNoSampleSection.csv"))
  }

  it should "throw an exception when a row in the sample section has fewer entries than expected" in {
    an[Exception] should be thrownBy SampleSheet(testDir.resolve("SampleSheetSampleMissingColumns.csv"))
  }

  it should "throw an exception when empty rows are mixed in with sample data rows" in {
    val sampleSheetLines = Io.readLines(testDir.resolve("SampleSheet.csv")).toSeq
    val sampleSheetPre   = sampleSheetLines.take(28) // gets us to Sample_Name_6
    val sampleSheetPost  = sampleSheetLines.drop(28) // Sample_Name_7 and later

    // Sample sheet with empty lines in the middle of sample data section
    val emptyLines     = sampleSheetPre.iterator ++ Seq("", "") ++ sampleSheetPost
    val emptyLineError = the [IllegalArgumentException] thrownBy SampleSheet(lines = emptyLines, None)
    emptyLineError.getMessage should startWith ("Found empty line on row")

    // Sample sheet with comma-only lines in the middle of sample data section
    val commaLines     = sampleSheetPre.iterator ++ Seq(",,,,,,,", ",,,,,,,,") ++ sampleSheetPost
    val commaLineError = the [IllegalArgumentException] thrownBy SampleSheet(lines = commaLines, lane = None)
    commaLineError.getMessage should startWith ("Found empty line on row")

    // Sample sheet with mix of empty and comma-only lines in the middle of sample data section
    val commaAndEmptyLines     = sampleSheetPre.iterator ++ Seq(",,,,,,,", ",,,,,,,,", "", "", ",,,,,,,,", "") ++ sampleSheetPost
    val commaAndEmptyLineError = the [IllegalArgumentException] thrownBy SampleSheet(lines = commaAndEmptyLines, lane = None)
    commaAndEmptyLineError.getMessage should startWith ("Found empty line on row")
  }

  it should "throw an exception if sample ids are not unique " in {
    // with no lane defined
    val sample = Sample(sampleOrdinal=0, sampleId="ID", sampleName="NAMEA", libraryId="LIBRARY")
    an[IllegalArgumentException] should be thrownBy new SampleSheet(Seq(sample, sample.copy(sampleName="NAMEB")))

    // within a lane
    val sampleForLane = Sample(sampleOrdinal=0, sampleId="ID", sampleName="NAMEA", libraryId="LIBRARY", lane=Some(1))
    an[IllegalArgumentException] should be thrownBy new SampleSheet(Seq(sampleForLane, sampleForLane.copy(sampleName="NAMEB")))

    // where one sample defines lane, another does not
    val sampleWithLane = Sample(sampleOrdinal=0, sampleId="ID", sampleName="NAMEA", libraryId="LIBRARY", lane=Some(1))
    val sampleNoLane   = sampleWithLane.copy(sampleName="NAMEB", lane=None)
    an[IllegalArgumentException] should be thrownBy new SampleSheet(Seq(sampleWithLane, sampleNoLane))
  }

  it should "throw an exception if the combination of sample name and library id are not unique" in {
    // Different sample ID, same sample name and library ID, lane is None
    val sample = Sample(sampleOrdinal=0, sampleId="ID1", sampleName="NAME", libraryId="LIBRARY")
    an[IllegalArgumentException] should be thrownBy new SampleSheet(Seq(sample, sample.copy(sampleId="ID2")))

    // Different sample ID, same sample name and library ID, lane is 1
    val sampleForLane = Sample(sampleOrdinal=0, sampleId="ID1", sampleName="NAME", libraryId="LIBRARY", lane=Some(1))
    an[IllegalArgumentException] should be thrownBy new SampleSheet(Seq(sampleForLane, sampleForLane.copy(sampleId="ID2")))
  }

  it should "not throw an exception if sample ids are unique and if the combination of sample name and library id are unique" in {
    // Different sample ID, different sample name, same library ID, no lane
    {
      val sample = Sample(sampleOrdinal = 0, sampleId = "IDA", sampleName = "NAMEA", libraryId = "LIBRARYA")
      new SampleSheet(Seq(sample, sample.copy(sampleId = "IDB", sampleName = "NAMEB"))).size shouldBe 2
    }

    // Different sample ID, different sample name, same library ID, same lane
    {
      val sample = Sample(sampleOrdinal = 0, sampleId = "IDA", sampleName = "NAMEA", libraryId = "LIBRARYA", lane = Some(1))
      new SampleSheet(Seq(sample, sample.copy(sampleId = "IDB", sampleName = "NAMEB"))).size shouldBe 2
    }

    // Different sample ID, same sample name, different library ID, same lane
    {
      val sample = Sample(sampleOrdinal = 0, sampleId = "IDA", sampleName = "NAMEA", libraryId = "LIBRARYA", lane = Some(1))
      new SampleSheet(Seq(sample, sample.copy(sampleId = "IDB", libraryId = "LIBRARYB"))).size shouldBe 2
    }

    // Different sample ID, same sample name, same library ID, different lane
    {
      val sample = Sample(sampleOrdinal = 0, sampleId = "IDA", sampleName = "NAMEA", libraryId = "LIBRARY", lane = Some(1))
      new SampleSheet(Seq(sample, sample.copy(sampleId = "IDB", lane = Some(2)))).size shouldBe 2
    }
  }

  it should "not throw an exception if samples are replicated across lanes" in {
    val sample1 = Sample(sampleOrdinal=0, sampleId="ID1", sampleName="NAME1", libraryId="LIBRARY1", lane=Some(1))
    val sample2 = Sample(sampleOrdinal=1, sampleId="ID2", sampleName="NAME2", libraryId="LIBRARY2", lane=Some(1))

    val samples = Seq(
      sample1,
      sample2,
      sample1.copy(lane=Some(2)),
      sample2.copy(lane=Some(2))
    )

    new SampleSheet(samples).size shouldBe 4
  }

  "SampleSheet.get" should "get the ith sample" in {
    val sampleSheet = SampleSheet(testDir.resolve("SampleSheet.csv"))
    val samples = sampleSheet.toList
    samples.zipWithIndex.foreach { case (sample, i) =>
        sampleSheet.get(i) shouldBe samples(i)
    }
  }

  "SampleSheet.cleanMiseqSampleId" should "clean [#] to ' ' and [_+ ] to '-'" in {
    SampleSheet.cleanMiseqSampleId("ABCDEFGHIJKLMNOP") shouldBe "ABCDEFGHIJKLMNOP"
    SampleSheet.cleanMiseqSampleId("#") shouldBe " "
    SampleSheet.cleanMiseqSampleId("_+ ") shouldBe "---"
    SampleSheet.cleanMiseqSampleId("ABCDE#_+ ") shouldBe "ABCDE ---"
  }

  "SampleSheet.SplitRegex" should "split by commas but commas in quotes" in {
    "A,B,C".split(",") should have size 3
    "A,B,C".split(SampleSheet.SplitRegex) should have size 3
    "A,\"B,C\",D".split(SampleSheet.SplitRegex) should have size 3
    "\"B,C\",D,E".split(SampleSheet.SplitRegex) should have size 3
    "A,B,\"C,D\"".split(SampleSheet.SplitRegex) should have size 3
  }
}

