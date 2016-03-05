package com.fulcrumgenomics.fastq

import com.fulcrumgenomics.testing.UnitSpec
import dagr.commons.io.Io

object FastqSourceTest {
  val someFastq =
    """
      |@Foo:279:000000000-ABCDE:1:1101:15033:1749 1:N:0:18
      |ATGACTGCGTATATATGACGTCGTGCTA
      |+
      |-,86,,;:C,<;-,86,,;:C,<;-,86
      |@Foo:279:000000000-ABCDE:1:1101:12430:1790 1:N:0:18
      |ATAGCTATGAGCTGCTGCTGA
      |+
      |8,,8,++++8788,,8,++++
      |@Foo:279:000000000-ABCDE:1:1101:9471:1291/1 1:N:0:18
      |CCCCCCCTCCTTGGGGGGGAGAGGG
      |+Foo:279:000000000-ABCDE:1:1101:9471:1291/1 1:N:0:18
      |AAA,,++78,66AAA,,++78,666
    """.stripMargin.trim().lines.toSeq.dropWhile(_.isEmpty)

  val someFastqRecords = Seq(
    FastqRecord("Foo:279:000000000-ABCDE:1:1101:15033:1749", "ATGACTGCGTATATATGACGTCGTGCTA", "-,86,,;:C,<;-,86,,;:C,<;-,86", Some("1:N:0:18"), None),
    FastqRecord("Foo:279:000000000-ABCDE:1:1101:12430:1790", "ATAGCTATGAGCTGCTGCTGA",        "8,,8,++++8788,,8,++++",        Some("1:N:0:18"), None),
    FastqRecord("Foo:279:000000000-ABCDE:1:1101:9471:1291",  "CCCCCCCTCCTTGGGGGGGAGAGGG",    "AAA,,++78,66AAA,,++78,666",    Some("1:N:0:18"), Some(1))
  )
}

/**
  * Tests for the FastqSource
  */
class FastqSourceTest extends UnitSpec {
  "FastqSource" should "read valid fastq from various kinds of input resources" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, FastqSourceTest.someFastq)

    Seq(FastqSource(fq), FastqSource(fq.toFile), FastqSource(Io.toInputStream(fq)), FastqSource(FastqSourceTest.someFastq)).foreach(source => {
      source.hasNext shouldBe true
      source.next shouldBe FastqSourceTest.someFastqRecords(0)
      source.hasNext shouldBe true
      source.next shouldBe FastqSourceTest.someFastqRecords(1)
      source.hasNext shouldBe true
      source.next shouldBe FastqSourceTest.someFastqRecords(2)
      source.hasNext shouldBe false
      an[NoSuchElementException] shouldBe thrownBy { source.next() }
    })
  }

  it should "throw an exception if a header line doesn't start with an @" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, Seq("foo", "ACGT", "+", "7777"))
    an[IllegalStateException] shouldBe thrownBy { FastqSource(fq) }
  }

  it should "throw an exception if a quality header line doesn't start with an +" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, Seq("@foo", "ACGT", "@", "7777"))
    an[IllegalStateException] shouldBe thrownBy { FastqSource(fq) }
  }

  it should "throw an exception if there are not n*4 lines" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, Seq("@foo", "ACGT", "+foo"))
    an[IllegalStateException] shouldBe thrownBy { FastqSource(fq) }
  }

  it should "throw an exception if seqs and quals are not the same length for a record" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, Seq("@foo", "ACGT", "+foo", "12345"))
    an[IllegalStateException] shouldBe thrownBy { FastqSource(fq) }
  }
}
