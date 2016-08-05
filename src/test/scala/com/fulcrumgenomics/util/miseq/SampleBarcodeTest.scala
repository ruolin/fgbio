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

package com.fulcrumgenomics.util.miseq

import com.fulcrumgenomics.testing.UnitSpec

/**
  * Tests for SampleBarcode.
  */
class SampleBarcodeTest extends UnitSpec {
  "SampleBarcode.hashCode" should "be the hash code of the concatenated barcode sequences" in {
    new SampleBarcode(Seq("GATTACA")).hashCode shouldBe "GATTACA".hashCode
    new SampleBarcode(Seq("A", "C", "G", "T")).hashCode shouldBe "A-C-G-T".hashCode
  }

  "SampleBarcode.equals" should "compare based on the concatenated barcode sequences" in {
    new SampleBarcode(Seq("GATTACA")).equals(2) shouldBe false
    new SampleBarcode(Seq("GATTACA")).equals(new SampleBarcode(Seq("GATTACA"))) shouldBe true
    new SampleBarcode(Seq("GATTACA")).equals(new SampleBarcode(Seq("GATTACT"))) shouldBe false
  }

  "SampleBarcode.toString" should "be the concatened barcode" in {
    new SampleBarcode(Seq("GATTACA")).toString shouldBe "GATTACA"
    new SampleBarcode(Seq("A", "C", "G", "T")).toString shouldBe "A-C-G-T"
  }
}
