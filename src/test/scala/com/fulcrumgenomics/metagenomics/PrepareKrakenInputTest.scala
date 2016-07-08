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

package com.fulcrumgenomics.metagenomics

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fastq.{FastqRecord, FastqWriter}
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.SAMUtils


/** Tests for PrepareKrakenInput. */
class PrepareKrakenInputTest extends UnitSpec {
  /** Creates a fastq record with arbitrary bases, and qualities described by the qs
    * parameter, which is interpreted as (qscore, count).
    */
  def fq(name: String, qs: (Int,Int)*): FastqRecord = {
    val quals = qs.flatMap { case (q, n) => List.fill(n)(SAMUtils.phredToFastq(q)) }.mkString
    val bases = quals.map(ch => 'A')
    FastqRecord(name, bases=bases, quals=quals)
  }

  case class FastaRecord(header: String, bases: String)

  /** Runs the command on either a single or pair of fastq records and returns the result. */
  def exec(r1: FastqRecord, r2: Option[FastqRecord], q:Int=20, len:Int=30): Option[FastaRecord] = {
    val fq1   = makeTempFile("prepare_kraken_input_test.", ".fq")
    val fq2   = makeTempFile("prepare_kraken_input_test.", ".fq")
    val fasta = makeTempFile("prepare_kraken_input_test.", ".fasta")

    FastqWriter(fq1).write(r1).close()
    r2.foreach(r => FastqWriter(fq2).write(r).close())
    new PrepareKrakenInput(read1=Seq(fq1), read2=if (r2.isEmpty) Seq() else Seq(fq2), output=fasta, minQuality=q, minReadLength=len, kmerSize=20).execute()

    val lines = Io.readLines(fasta).toSeq
    fq1.toFile.delete()
    fq2.toFile.delete()
    fasta.toFile.delete()

    lines match {
      case Seq(header, bases) => Some(FastaRecord(header, bases))
      case Seq()              => None
      case _                  => unreachable("should not have more than two lines of output!")
    }
  }


  /** Tests for the method that just trims and masks a single fastq record. */
  {
    "PrepareKrakenInput.trimAndMask" should "not trim or mask anything in" in {
      val prepper = new PrepareKrakenInput(read1=Seq(), read2=Seq(), minQuality=20, minReadLength=10)
      val bases = prepper.trimAndMask(fq("foo", (30, 20)))
      bases should have length 20
      bases should not contain 'N'
    }

    it should "trim but not mask" in {
      val prepper = new PrepareKrakenInput(read1=Seq(), read2=Seq(), minQuality=20, minReadLength=10)
      val bases = prepper.trimAndMask(fq("foo", (30, 20), (10, 20)))
      bases should have length 20
      bases should not contain 'N'
    }

    it should "mask one base but not trim" in {
      val prepper = new PrepareKrakenInput(read1=Seq(), read2=Seq(), minQuality=20, minReadLength=10)
      val bases = prepper.trimAndMask(fq("foo", (30, 5), (10, 1), (25, 5), (2, 2), (30,7)))
      bases should have length 20
      bases.count(_ == 'N') shouldBe 3
    }

    it should "trim and mask one base" in {
      val prepper = new PrepareKrakenInput(read1=Seq(), read2=Seq(), minQuality=20, minReadLength=10)
      val bases = prepper.trimAndMask(fq("foo", (5, 1), (30, 19), (10, 20)))
      bases should have length 20
      bases(0) shouldBe 'N'
      bases.drop(1) should not contain 'N'
    }
  }


  /** Tests for the command line program. */
  {
    "PrepareKrakenInput" should "output a single read without any modifications" in {
      val rec = exec(r1=fq("t1", (30, 100)), r2=None)
      rec shouldBe 'defined
      rec.get.header shouldBe ">t1"
      rec.get.bases should have length 100
      rec.get.bases should not contain 'N'
    }

    "it" should "output a pair of reads concatenated without any other modifications" in {
      val rec = exec(r1=fq("t1", (30, 50)), r2=Some(fq("t1", (30, 50))))
      rec shouldBe 'defined
      rec.get.header shouldBe ">t1"
      rec.get.bases should have length 120
    }

    "it" should "output nothing if the input reads are too short" in {
      exec(r1=fq("t1", (30, 10)), r2=Some(fq("t1", (30, 10)))) shouldBe None // Both too short
      exec(r1=fq("t1", (30, 100)), r2=Some(fq("t1", (30, 10)))) shouldBe None // R1 too short
      exec(r1=fq("t1", (30, 10)), r2=Some(fq("t1", (30, 100)))) shouldBe None // R2 too short
    }

    "it" should "output nothing if the input reads are too short post-trimming" in {
      exec(r1=fq("t1", (10, 100)), r2=Some(fq("t1", (10, 100)))) shouldBe None
    }

    "it" should "fail if the fastqs are out of sync" in {
      an[Exception] should be thrownBy exec(r1=fq("t1", (30, 10)), r2=Some(fq("t2", (30, 10))))
    }
  }
}
