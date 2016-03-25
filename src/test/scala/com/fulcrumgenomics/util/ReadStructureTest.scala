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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.UnitSpec

import scala.util.Try

class ReadStructureTest extends UnitSpec {

  private def compareReadStructures(actual: ReadStructure, expected: Seq[ReadSegment]): Unit = {
    // make sure the segements match
    actual.toSeq shouldBe expected
    // make sure the string representations are the same
    actual.toString shouldBe new ReadStructure(expected).toString
  }

  "ReadStructure" should "be built from a string" in {
    compareReadStructures(ReadStructure("1T"), Seq(Template(0, 1)))
    compareReadStructures(ReadStructure("1B"), Seq(SampleBarcode(0, 1)))
    compareReadStructures(ReadStructure("1M"), Seq(MolecularBarcode(0, 1)))
    compareReadStructures(ReadStructure("1S"), Seq(Skip(0, 1)))
    compareReadStructures(ReadStructure("101T"), Seq(Template(0, 101)))
    compareReadStructures(ReadStructure("5B101T"), Seq(SampleBarcode(0, 5), Template(5, 101)))
    compareReadStructures(ReadStructure("123456789T"), Seq(Template(0, 123456789)))
    compareReadStructures(ReadStructure("10T10B10B10S10M"), Seq(Template(0, 10), SampleBarcode(10, 10), SampleBarcode(20, 10), Skip(30, 10), MolecularBarcode(40, 10)))
  }

  it should "not be built from invalid structures" in {
    Try { ReadStructure("0T")        }.failed.get.getMessage should include ("[0]T")
    Try { ReadStructure("9R")        }.failed.get.getMessage should include ("9[R]")
    Try { ReadStructure("T")         }.failed.get.getMessage should include ("[T]")
    Try { ReadStructure("23TT")      }.failed.get.getMessage should include ("23T[T]")
    Try { ReadStructure("23T2")      }.failed.get.getMessage should include ("23T[2]")
    Try { ReadStructure("23T2TT23T") }.failed.get.getMessage should include ("23T2T[T]23T")
  }

  "ReadStructure.totalBases" should "return the total bases described by the read structure" in {
    ReadStructure("1T").totalBases shouldBe 1
    ReadStructure("1B").totalBases shouldBe 1
    ReadStructure("1M").totalBases shouldBe 1
    ReadStructure("1S").totalBases shouldBe 1
    ReadStructure("101T").totalBases shouldBe 101
    ReadStructure("123456789T").totalBases shouldBe 123456789
    ReadStructure("10T10B10B10S10M").totalBases shouldBe 50
  }

  "ReadStructure.extract" should "get extract the bases for each segment" in {
    val rs = ReadStructure("2T2B2M2S")
    rs.structureRead("AACCGGTT").foreach {
      case (bases: String, segment: ReadSegment) =>
        segment match {
          case t: Template         => bases shouldBe "AA"
          case s: SampleBarcode    => bases shouldBe "CC"
          case m: MolecularBarcode => bases shouldBe "GG"
          case x: Skip             => bases shouldBe "TT"
        }
    }
    an[IllegalStateException] should be thrownBy rs.structureRead("AAAAAAA")
  }
}