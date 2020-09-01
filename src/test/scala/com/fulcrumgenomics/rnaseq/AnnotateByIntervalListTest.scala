package com.fulcrumgenomics.rnaseq

import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
import com.fulcrumgenomics.testing.{ErrorLogLevel, ReferenceSetBuilder, SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.{IntervalListWriter, Metric}
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.Interval
import org.scalatest.Matchers.convertToAnyShouldWrapper
import org.scalatest.OptionValues

class AnnotateByIntervalListTest extends UnitSpec with ErrorLogLevel with OptionValues {
  private val (ref, refPath) = {
    val builder = new ReferenceSetBuilder
    builder.add("chr1").add("AAAAAAAAAA", 50) // 500 bases
    builder.add("chr2").add("CCCCCCCCCC", 50) // 500 bases

    val refPath = builder.toTempFile()
    val ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(refPath)
    (ref, refPath)
  }

  val dict: SAMSequenceDictionary = ref.getSequenceDictionary

  "AnnotateByIntervalList" should "associate template bases with interval if read one overlaps a single interval" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 150, false, "interval1")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 50
  }

  "AnnotateByIntervalList" should "associate template bases with interval if read two overlaps a single interval" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 200, 350, false, "interval1")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 50
  }

  "AnnotateByIntervalList" should "associate template bases with interval if both reads overlap a single interval" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 350, false, "interval1")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 100
  }

  "AnnotateByIntervalList" should "associate partial template bases with interval when part of a read overlaps" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 76, 275, false, "interval1")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    // 25 from read one 25 from read two
    metrics.filter(_.name == "interval1").head.bases shouldBe 50
  }

  "AnnotateByIntervalList" should "associate bases with abutting intervals" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 76, false, "interval1")
    writer += new Interval("chr1", 77, 100, false, "interval2")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    // 25 from read one 25 from read two
    metrics.filter(_.name == "interval1").head.bases shouldBe 50
  }
}
