/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.testing.UnitSpec
import dagr.commons.io.PathUtil
import dagr.commons.CommonsDef._
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.JavaConversions._

/**
  * Tests for HapCutToVcf.
  */
class HapCutToVcfTest extends UnitSpec {

  private val dir         = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/vcf")
  private val originalVcf = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.vcf")
  private val hapCut      = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut")
  private val output      = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut.vcf")
  private val outputGatk  = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut.gatk.vcf")

  private def countVcfRecords(vcf: PathToVcf): Int = {
    val vcfReader = new VCFFileReader(vcf.toFile, false)
    yieldAndThen(vcfReader.iterator().length)(vcfReader.close())
  }

  private def compareVcfs(newVcf: PathToVcf, originalVcf: PathToVcf): Unit = {
    val newVcfReader      = new VCFFileReader(newVcf.toFile, false)
    val originalVcfReader = new VCFFileReader(originalVcf.toFile, false)

    for (newVariantCtx <- newVcfReader) {
      originalVcfReader.exists { originalVariantCtx =>
        originalVariantCtx.getContig == newVariantCtx.getContig &&
        originalVariantCtx.getStart  == newVariantCtx.getStart &&
        originalVariantCtx.getEnd    == newVariantCtx.getEnd
      } shouldBe true
    }
  }

  private def getNumPhasedFromVcf(path: PathToVcf, gatkPhasingFormat: Boolean): Int = {
    val vcfReader = new VCFFileReader(path.toFile, false)
    val numPhased = vcfReader.iterator().count { ctx =>
      if (gatkPhasingFormat) ctx.isNotFiltered // how many are marked as passed filter
      else ctx.getGenotypes.exists(_.isPhased) // how many are marked as phased
    }
    vcfReader.close()
    numPhased
  }

  "HapCutReader" should "read in a HAPCUT file" in {
    val reader = new HapCutReader(hapCut)
    val allCalls = reader.toSeq
    val calls = allCalls.flatten
    allCalls.length shouldBe 342 // 342 total variants
    calls.length shouldBe 237 // 237 phased variants
    calls.map(_.phaseSet).distinct.length shouldBe 8 // eight phased blocks
  }

  "HapCutToVcf" should "convert a HAPCUT file to VCF in both GATK and VCF-spec phasing format" in {
    Stream(true, false).foreach { gatkPhasingFormat =>
      val expectedOutput = if (gatkPhasingFormat) outputGatk else output
      val out = makeTempFile("hap_cut_to_vcf.hapcut", ".vcf")

      new HapCutToVcf(
        vcf               = originalVcf,
        input             = hapCut,
        output            = out,
        gatkPhasingFormat = gatkPhasingFormat
      ).execute()

      // check that we have the same # of records in the output as the input
      countVcfRecords(out) shouldBe countVcfRecords(originalVcf)

      // check that all records in the output are found in the input
      compareVcfs(out, originalVcf)

      // get the # of phased variants from the output
      val numPhasedFromOut = getNumPhasedFromVcf(out, gatkPhasingFormat)

      // check that the # of variants phased in the output agrees with the # of phased calls produced by HapCut
      val hapCutReader = new HapCutReader(hapCut)
      val numPhasedFromHapCut = hapCutReader.flatten.length
      numPhasedFromOut shouldBe numPhasedFromHapCut

      // check that the # of variants phased in the output agrees with the # of phased calls in the expected output
      numPhasedFromOut shouldBe getNumPhasedFromVcf(expectedOutput, gatkPhasingFormat)
    }
  }
}
