package com.fulcrumgenomics.rnaseq

import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
import com.fulcrumgenomics.testing.{ErrorLogLevel, ReferenceSetBuilder, SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.{IntervalListWriter, Metric}
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.Interval
import org.scalatest.OptionValues

class AnnotateByIntervalListTest extends UnitSpec with ErrorLogLevel with OptionValues {
  private val ref = {
    val builder = new ReferenceSetBuilder
    builder.add("chr1").add("AAAAAAAAAA", 50) // 500 bases
    builder.add("chr2").add("CCCCCCCCCC", 50) // 500 bases

    val ref     = ReferenceSequenceFileFactory.getReferenceSequenceFile(builder.toTempFile())
    ref
  }

  val dict: SAMSequenceDictionary = ref.getSequenceDictionary

  "AnnotateByIntervalList" should "associate template bases if the template is enclosed in an interval" in {
    val builder       = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath  = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput     = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer        = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 350, false, "interval1")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 150
  }

  "AnnotateByIntervalList" should "not associate template bases if the template is not fully enclosed in an interval" in {
    val builder       = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath  = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput     = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer        = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 200, false, "interval1")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 0
  }

  "AnnotateByIntervalList" should "associate bases with only enclosing intervals" in {
    val builder       = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath  = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput     = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer        = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 500, false, "interval1")
    writer += new Interval("chr1", 152, 300, false, "interval2")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput, nameRegex = Some("(^interval.).*")).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 150
    metrics.filter(_.name == "interval2").head.bases shouldBe 0
  }

  "AnnotateByIntervalList" should "associate bases with all enclosing intervals (duplicate counts)" in {
    val builder       = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath  = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput     = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer        = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 500, false, "interval1")
    writer += new Interval("chr1", 10, 490, false, "interval2")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput, nameRegex = Some("(^interval.).*")).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 150
    metrics.filter(_.name == "interval1").head.unique_bases shouldBe 0
    metrics.filter(_.name == "interval2").head.bases shouldBe 150
    metrics.filter(_.name == "interval2").head.unique_bases shouldBe 0
  }

  "AnnotateByIntervalList" should "associate bases with the smallest interval (most specific counts)" in {
    val builder       = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath  = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput     = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer        = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 500, false, "interval1")
    writer += new Interval("chr1", 10, 490, false, "interval2")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput, nameRegex = Some("(^interval.).*"),
      overlapAssociationStrategy = OverlapAssociationStrategy.MostSpecific).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 0
    metrics.filter(_.name == "interval1").head.unique_bases shouldBe 0
    metrics.filter(_.name == "interval2").head.bases shouldBe 150
    metrics.filter(_.name == "interval2").head.unique_bases shouldBe 150
  }

  "AnnotateByIntervalList" should "associate bases with the largest interval (least specific counts)" in {
    val builder       = new SamBuilder(readLength=50, sd=Some(dict.fromSam))
    val intervalPath  = makeTempFile("annotatebyinterval", "interval_list")
    val tmpOutput     = makeTempFile("annotatebyinterval_metrics", "txt")
    val writer        = IntervalListWriter(intervalPath, dict.fromSam)

    writer += new Interval("chr1", 1, 500, false, "interval1")
    writer += new Interval("chr1", 10, 490, false, "interval2")
    writer.close()

    builder.addPair("single", start1 = 50, start2 = 250)
    new AnnotateByIntervalList(builder.toTempFile(), intervalPath, tmpOutput, nameRegex = Some("(^interval.).*"),
      overlapAssociationStrategy = OverlapAssociationStrategy.LeastSpecific).execute()

    val metrics = Metric.read[CoverageByIntervalMetrics](tmpOutput)
    metrics.filter(_.name == "interval1").head.bases shouldBe 150
    metrics.filter(_.name == "interval2").head.bases shouldBe 0
  }
}
