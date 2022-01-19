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

package com.fulcrumgenomics.bam

import java.nio.file.Paths
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.{Amplicon, Io, Metric}
import htsjdk.samtools.util.{CoordMath, SequenceUtil}
import htsjdk.samtools.{CigarOperator => Op}

class TrimPrimersTest extends UnitSpec {
  val amplicons: Seq[Amplicon] = Seq(
    Amplicon("chr1", 100, 119, 200, 219),
    Amplicon("chr1", 200, 219, 300, 320),
    Amplicon("chr1", 300, 319, 400, 421),
    Amplicon("chr1", 400, 419, 500, 522)
  )

  val refFasta = Paths.get("src/test/resources/com/fulcrumgenomics/bam/trim_primers_test.fa")

  val maxPrimerLength = CoordMath.getLength(500, 522)
  val readLength = 50

  def primers: FilePath = {
    val tmp = makeTempFile("primers.", ".txt")
    Metric.write(path=tmp, amplicons)
    tmp
  }

  "TrimPrimers" should "trim reads that match primer locations" in {
    val builder = new SamBuilder(readLength=readLength, sort=Some(SamOrder.Coordinate))
    builder.addPair("q1", start1=100, start2=CoordMath.getStart(219, readLength))
    builder.addPair("q2", start1=200, start2=CoordMath.getStart(320, readLength))
    builder.addPair("q3", start1=300, start2=CoordMath.getStart(421, readLength))
    builder.addPair("q4", start1=400, start2=CoordMath.getStart(522, readLength))
    val bam = builder.toTempFile()
    val newBam = makeTempFile("trimmed.", ".bam")
    new TrimPrimers(input=bam, output=newBam, primers=primers, hardClip=false).execute()

    val reads = readBamRecs(newBam)
    reads should have size 8
    reads.foreach { rec =>
      if (rec.negativeStrand) rec.cigar.elems.last.operator shouldBe Op.SOFT_CLIP
      else                    rec.cigar.elems.head.operator shouldBe Op.SOFT_CLIP
    }
  }

  it should "trim reads that are off from primer locations by a little bit" in {
    val builder = new SamBuilder(readLength=readLength, sort=Some(SamOrder.Coordinate))
    builder.addPair("q1", start1=100+2, start2=CoordMath.getStart(219, readLength)-1)
    builder.addPair("q2", start1=200+1, start2=CoordMath.getStart(320, readLength)+2)
    builder.addPair("q3", start1=300+2, start2=CoordMath.getStart(421, readLength)+1)
    builder.addPair("q4", start1=400-1, start2=CoordMath.getStart(522, readLength)-2)
    val bam = builder.toTempFile()
    val newBam = makeTempFile("trimmed.", ".bam")
    new TrimPrimers(input=bam, output=newBam, primers=primers, hardClip=false).execute()

    val reads = readBamRecs(newBam)
    reads should have size 8
    reads.foreach { rec =>
      if (rec.negativeStrand) rec.cigar.elems.last.operator shouldBe Op.SOFT_CLIP
      else                    rec.cigar.elems.head.operator shouldBe Op.SOFT_CLIP
    }
  }

  it should "trim FR reads regardless of which is first/second of pair" in {
    val builder = new SamBuilder(readLength=readLength, sort=Some(SamOrder.Coordinate))
    builder.addPair("q1", start2=100, strand2=Plus, start1=CoordMath.getStart(219, readLength), strand1=Minus)
    builder.addPair("q2", start2=200, strand2=Plus, start1=CoordMath.getStart(320, readLength), strand1=Minus)
    builder.addPair("q3", start2=300, strand2=Plus, start1=CoordMath.getStart(421, readLength), strand1=Minus)
    builder.addPair("q4", start2=400, strand2=Plus, start1=CoordMath.getStart(522, readLength), strand1=Minus)
    val bam = builder.toTempFile()
    val newBam = makeTempFile("trimmed.", ".bam")
    new TrimPrimers(input=bam, output=newBam, primers=primers, hardClip=false).execute()

    val reads = readBamRecs(newBam)
    reads should have size 8
    reads.foreach { rec =>
      if (rec.negativeStrand) rec.cigar.elems.last.operator shouldBe Op.SOFT_CLIP
      else                    rec.cigar.elems.head.operator shouldBe Op.SOFT_CLIP
    }
  }

  it should "trim anything that's not an on-amplicon FR pair by the longest primer length" in {
    val builder = new SamBuilder(readLength=readLength, sort=Some(SamOrder.Coordinate))
    builder.addPair("q1", start1=100, start2=CoordMath.getStart(219, readLength)+20)            // Too far from primer sites
    builder.addPair("q2", start1=200, start2=CoordMath.getStart(320, readLength), strand2=Plus) // FF pair
    builder.addPair("q3", start1=300, start2=300, unmapped2=true)                               // Unmapped mate
    builder.addFrag("q4", start=CoordMath.getStart(219, readLength), strand=Minus)              // Fragment read
    val bam = builder.toTempFile()
    val newBam = makeTempFile("trimmed.", ".bam")
    new TrimPrimers(input=bam, output=newBam, primers=primers, hardClip=false).execute()

    val reads = readBamRecs(newBam)
    reads should have size 7
    reads.filter(r => r.mapped).foreach { rec =>
      val elem = if (rec.negativeStrand) rec.cigar.elems.last else rec.cigar.elems.head
      elem.operator shouldBe Op.SOFT_CLIP
      elem.length   shouldBe maxPrimerLength
    }
  }

  it should "trim back reads that are more than fully overlapped after clipping" in {
    val builder = new SamBuilder(readLength=120, sort=Some(SamOrder.Coordinate))
    builder.addPair("q1", start1=100, start2=100)
    val bam = builder.toTempFile()
    val newBam = makeTempFile("trimmed.", ".bam")
    new TrimPrimers(input=bam, output=newBam, primers=primers, hardClip=false).execute()

    val reads = readBamRecs(newBam)
    reads should have size 2
    reads.foreach ( rec => rec.cigar.toString shouldBe "20S80M20S" )
  }

  it should "recalculate NM/UQ/MD if a reference is given" in {
    import SamOrder._
    val zeroErrors   = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    val oneErrors    = "AAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    val twoErrors    = "AAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAA"
    val threeErrors  = "AAAAAGAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAGAAAAA"
    def rc(s: String) = SequenceUtil.reverseComplement(s)

    val orders = Seq(Some(Queryname), Some(Coordinate), None)
    for (inOrder <- orders; outOrder <- orders) {
      val builder = new SamBuilder(readLength=readLength, sort=inOrder)
      builder.addPair("q1", start1=100, start2=CoordMath.getStart(219, readLength), bases1=zeroErrors,  bases2=zeroErrors.reverse )
      builder.addPair("q2", start1=200, start2=CoordMath.getStart(320, readLength), bases1=oneErrors,   bases2=oneErrors.reverse  )
      builder.addPair("q3", start1=300, start2=CoordMath.getStart(421, readLength), bases1=twoErrors,   bases2=twoErrors.reverse  )
      builder.addPair("q4", start1=400, start2=CoordMath.getStart(522, readLength), bases1=threeErrors, bases2=threeErrors.reverse)

      val bam = builder.toTempFile()
      val newBam = makeTempFile("trimmed.", ".bam")
      new TrimPrimers(input=bam, output=newBam, primers=primers, hardClip=false, ref=Some(refFasta), sortOrder=outOrder).execute()

      val reads = readBamRecs(newBam)
      reads should have size 8
      reads.foreach { rec =>
        rec.get[Int]("NM")    shouldBe defined
        rec.get[String]("MD") shouldBe defined
        rec.get[Int]("UQ")    shouldBe defined

        val expectedNm = (rec.negativeStrand, rec.name) match {
          case (false, "q1") => 0
          case (true , "q1") => 0
          case (false, "q2") => 0
          case (true , "q2") => 0
          case (false, "q3") => 1
          case (true , "q3") => 1
          case (false, "q4") => 2
          case (true , "q4") => 2
          case (_,        _) => unreachable()
        }

        rec[Int]("NM") shouldBe expectedNm
      }
    }
  }

  it should "fail if no primers are given in the primer file" in {
    val primerFile = makeTempFile("primers.", ".txt")
    val bam = new SamBuilder().toTempFile()
    Metric.write[Amplicon](primerFile)
    an[Exception] shouldBe thrownBy {
      new TrimPrimers(input=bam, output=bam, primers=primerFile, hardClip=false).execute()
    }
  }

  it should "fail if the headers are missing from the primer file" in {
    val primerFile = makeTempFile("primers.", ".txt")
    val bam = new SamBuilder().toTempFile()
    Io.writeLines(primerFile, Seq(Seq("chr1", 1000, 1020, 1200, 1220).mkString("\t")))
    an[Exception] shouldBe thrownBy {
      new TrimPrimers(input=bam, output=bam, primers=primerFile, hardClip=false).execute()
    }
  }

  it should "trim only R1s when --first-of-pair is given" in {
    val amplicons: Seq[Amplicon] = Seq(
      Amplicon("chr1", 100, 119, -1, -1),
      Amplicon("chr1", -1, -1, 300, 320),
      Amplicon("chr1", -1, -1, 601, 619),
    )
    val primers: FilePath = {
      val tmp = makeTempFile("primers.", ".txt")
      Metric.write(path=tmp, amplicons)
      tmp
    }

    val builder = new SamBuilder(readLength=readLength, sort=Some(SamOrder.Coordinate))

    // amplicon 1 - R1 should be trimmed (matches left primer)
    builder.addPair("q1", start1=100, start2=CoordMath.getStart(219, readLength), strand1=Plus, strand2=Minus)
    // amplicon 1 - R1 should be trimmed the maximum, since it matches the right primer (0-length)
    builder.addPair("q2", start2=100, start1=CoordMath.getStart(200, readLength), strand1=Minus, strand2=Plus)
    // amplicon 1 - R1 should be trimmed (matches left primer), since R2's coordinate is not considered
    builder.addPair("q3", start1=100, start2=500, strand1=Plus, strand2=Minus)

    // amplicon 2 - R1 should be trimmed (matches right primer)
    builder.addPair("q4", start1=CoordMath.getStart(619, readLength), start2=500, strand1=Minus, strand2=Plus)
    // amplicon 2 - R1 should be trimmed the maximum, since it matches the left primer (0-length)
    builder.addPair("q5", start2=CoordMath.getStart(619, readLength), start1=500, strand1=Plus, strand2=Minus)
    // amplicon 2 - R1 should be trimmed (matches right primer), since R2's coordinate is not considered
    builder.addPair("q6", start1=CoordMath.getStart(619, readLength), start2=400, strand1=Minus, strand2=Plus)

    val bam = builder.toTempFile()
    val newBam = makeTempFile("trimmed.", ".bam")
    new TrimPrimers(input=bam, output=newBam, primers=primers, hardClip=false, firstOfPair=true).execute()

    val reads = readBamRecs(newBam).sortBy(rec => (rec.name, rec.secondOfPair))
    reads should have size 12

    def validate(rec: SamRecord, name: String, firstOfPair: Boolean, fivePrimeSoftClipLength: Int = 0): Unit = {
      rec.name shouldBe name
      rec.firstOfPair shouldBe firstOfPair
      val elem = if (rec.negativeStrand) rec.cigar.elems.last else rec.cigar.elems.head
      if (rec.firstOfPair) {
        elem.operator shouldBe Op.SOFT_CLIP
        elem.length shouldBe fivePrimeSoftClipLength
      }
      else {
        elem.operator should not be Op.SOFT_CLIP
      }
    }

    validate(rec=reads(0),  name="q1", firstOfPair=true, fivePrimeSoftClipLength=20) // amplicon 1
    validate(rec=reads(1),  name="q1", firstOfPair=false)
    validate(rec=reads(2),  name="q2", firstOfPair=true, fivePrimeSoftClipLength=21) // maximum
    validate(rec=reads(3),  name="q2", firstOfPair=false)
    validate(rec=reads(4),  name="q3", firstOfPair=true, fivePrimeSoftClipLength=20) // amplicon 1
    validate(rec=reads(5),  name="q3", firstOfPair=false)
    validate(rec=reads(6),  name="q4", firstOfPair=true, fivePrimeSoftClipLength=19) // amplicon 2
    validate(rec=reads(7),  name="q4", firstOfPair=false)
    validate(rec=reads(8),  name="q5", firstOfPair=true, fivePrimeSoftClipLength=21) // maximum
    validate(rec=reads(9),  name="q5", firstOfPair=false)
    validate(rec=reads(10), name="q6", firstOfPair=true, fivePrimeSoftClipLength=19) // amplicon 2
    validate(rec=reads(11), name="q6", firstOfPair=false)
  }
}
