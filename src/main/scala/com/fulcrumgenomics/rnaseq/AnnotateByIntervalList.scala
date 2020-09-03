package com.fulcrumgenomics.rnaseq

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, IteratorToJavaCollectionsAdapter, javaIterableToIterator}
import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{FilePath, PathToBam, PathToIntervals}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
import com.fulcrumgenomics.rnaseq.OverlapAssociationStrategy.{Duplicate, LeastSpecific, MostSpecific}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Metric, ProgressLogger}
import enumeratum.EnumEntry
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.util.{CoordMath, Interval, IntervalList, OverlapDetector}

@clp(group = ClpGroups.Fastq, description =
  """
    |Collects base and read coverage statistics across intervals by name.
    |A regular expression may used to classify the name of the intervals to be counted and how they should be grouped.  Otherwise, the name is used as-is.
    |
    |For example for the following interval list:
    |
    |```
    |chr1    981971  982021  +       utr
    |chr1    982022  982065  +       intron
    |chr1    982066  982117  -       CDS;1
    |chr1    998963  999059  +       utr
    |chr1    999060  999432  -       CDS;2
    |chr1    999433  999526  +       intron
    |chr1    999527  999613  -       CDS;3
    |chr1    999614  999692  +       intron
    |chr1    999693  999973  -       CDS;4
    |```
    |
    |Using the the regex `([^;]*).*` would produce 3 separate interval counts CDS, utr and intron.
    |
    |Template bases (paired-end) or read bases (single end) are only counted if they fully enclosed within an interval.
    |By default if a read is enclosed by multiple intervals it will be counted for both intervals. This behavior
    |can be changed using a different `--overlap-association-strategy`.
  """)
class AnnotateByIntervalList
(
  @arg(flag = 'i', doc = "Path to a mapped BAM file.") val input: PathToBam,
  @arg(flag = 'l', doc = "Path to intervals") val intervals: PathToIntervals,
  @arg(flag = 'o', doc = "Path to the output metrics file") val output: FilePath,
  @arg(flag = 'R', doc = "The regular expression to use to classify the name of the interval") val nameRegex: Option[String] = None,
  @arg(flag = 'u', doc = "The name of the intervals not in the given set of intervals") val unmatched: String = "unmatched",
  @arg(flag = 'm', doc = "The minimum mapping quality for a read to be included") val minMappingQuality: Int = 20,
  @arg(flag = 'q', doc = "The minimum base quality for a base to be included") val minBaseQuality: Int = 0,
  @arg(flag = 'a', doc = "The strategy to determine which interval overlapping templates should be associated with")
  val overlapAssociationStrategy: OverlapAssociationStrategy = OverlapAssociationStrategy.Duplicate
) extends FgBioTool with LazyLogging {

  Io.assertReadable(this.input)
  Io.assertReadable(this.intervals)
  Io.assertCanWriteFile(this.output)

  override def execute(): Unit = {
    val progress = ProgressLogger(logger, verb = "Processed", noun = "loci", unit = 1000000)
    val reader = SamSource(this.input)
    val dict = reader.dict

    val (baseCounter, uniqueBaseCounter, templateCounter, uniqueTemplateCounter, overlapDetector) = {
      val _baseCounter            = new SimpleCounter[String]()
      val _uniqueBaseCounter      = new SimpleCounter[String]()
      val _templateCounter        = new SimpleCounter[String]()
      val _uniqueTemplateCounter  = new SimpleCounter[String]()
      val intervals               = buildIntervals(dict.asSam)
      intervals.foreach { interval => _baseCounter.count(interval.getName, 0) }

      val _overlapDetector = OverlapDetector.create[Interval](intervals.iterator.toJavaList)
      (_baseCounter, _uniqueBaseCounter, _templateCounter, _uniqueTemplateCounter, _overlapDetector)
    }

    // Templates are assigned to an interval together.
    val iterator = reader.iterator
    val records = iterator.filterNot(record => record.duplicate
      || record.secondary
      || record.supplementary
      || !record.pf
      || (record.paired && !record.isFrPair)
      || record.mapq < minMappingQuality)

    // First we filter to ensure that we only pick one read in a template and then iterate over those reads to create
    // an interval for each template. For single ended reads the interval is just the read.
    val templateIntervals = records.filter { record =>
      if (record.mateMapped) {
        val equalStartTieBreaker = {
          if (Math.min(record.start, record.end) == Math.min(record.mateStart, record.mateEnd.getOrElse(Int.MaxValue))) {
            record.firstOfPair
          } else {
            true
          }
        }
        // Take only the first read by coordinate order
        Math.min(record.start, record.end) < Math.min(record.mateStart, record.mateEnd.getOrElse(Int.MaxValue)) &&
        // if mate start == start pick first of pair
        equalStartTieBreaker
      } else {
        //If we are single ended then just take the record
        true
      }
    }).map(record =>
      // Interval should be end of read one to the beginning of read two
      new Interval(
        record.refName,
        Math.max(record.start, record.mateStart) + 1,
        Math.max(record.mateStart, record.mateEnd.getOrElse(exception) - 1
      )
    )

    var numLociMultiMatched: Long = 0

    templateIntervals.foreach { templateInterval =>
      // Overlaps should only be flagged when they enclose the template bases.
      val overlaps = overlapDetector.getOverlaps(templateInterval).toSeq
        .filter(interval => CoordMath.encloses(interval.getStart, interval.getEnd,
          templateInterval.getStart, templateInterval.getEnd))

      if (overlaps.isEmpty) {
        baseCounter.count(this.unmatched, templateInterval.getLengthOnReference)
        templateCounter.count(this.unmatched, 1)
      } else {
        val names: Seq[String] = overlaps.map(_.getName).distinct

        def countForNames(namesToCount: Seq[String]): Unit = {
          namesToCount.foreach { name =>
            val bases = overlaps.filter(_.getName == name).map(_.getIntersectionLength(templateInterval)).sum
            baseCounter.count(name, bases)
            templateCounter.count(name, 1)
            // If we only have one enclosing interval count unique bases and templates as well.
            if (namesToCount.length == 1) {
              uniqueBaseCounter.count(name, bases)
              uniqueTemplateCounter.count(name, 1)
            }
          }
        }

        if (names.length > 1) {
          // If there are multiple enclosing intervals by name use an associate strategy.
          overlapAssociationStrategy match {
            case Duplicate      => countForNames(names)
            case MostSpecific   => countForNames(Seq(overlaps.minBy(_.length).getName))
            case LeastSpecific  => countForNames(Seq(overlaps.maxBy(_.length).getName))
          }
          numLociMultiMatched += 1
        } else {
          // When there is only one interval matched we associate all intersected bases with the single interval
          countForNames(names)
        }
      }
      progress.record(templateInterval.getContig, templateInterval.getStart)
    }
    progress.logLast()
    reader.safelyClose()

    logger.info(f"Found $numLociMultiMatched%,d templates with multiple classifications.")

    // Get the total # of bases covered by each set of intervals
    val lengthByIntervalName = new SimpleCounter[String]()
    overlapDetector.getAll.toSeq.groupBy(_.getName).foreach { case (name, intvs) =>
      val intervalList = new IntervalList(dict.asSam)
      intervalList.addall(intvs.iterator.toJavaList)

      val length = intervalList.uniqued().map(_.length()).sum
      lengthByIntervalName.count(name, length)
    }
    lengthByIntervalName.count(this.unmatched, dict.length - lengthByIntervalName.total)

    val totalBases  = baseCounter.total.toDouble
    val metrics     = baseCounter.map { case (name, bases) =>
      // The total # of genomic bases covered by the intervals with the given name
      val span = lengthByIntervalName.countOf(name)

      CoverageByIntervalMetrics(
        name              = name,
        bases             = bases,
        unique_bases      = uniqueBaseCounter.countOf(name),
        span              = span,
        frac_genome       = span / dict.length.toDouble,
        frac_total_bases  = if (totalBases > 0) bases / totalBases else 0f,
        mean_coverage     = if (span > 0) bases / span.toDouble else 0f,
        templates         = templateCounter.countOf(name),
        unique_templates  = uniqueTemplateCounter.countOf(name)
      )
    }.toSeq.sortBy(_.name)
    Metric.write(this.output, metrics)
  }

  private def buildIntervals(dict: SAMSequenceDictionary): List[Interval] = {
    val inputIntervalList = IntervalList.fromPath(this.intervals)
    this.nameRegex.map(_.r) match {
      case None        => inputIntervalList.uniqued(true).toList
      case Some(regex) =>
        val intervals = inputIntervalList.map { i =>
          val name = i.getName match {
            case regex(n) => n
            case _        => throw new IllegalArgumentException(s"Could not match name '${i.getName}' with regex '${this.nameRegex}'")
          }
          new Interval(i.getContig, i.getStart, i.getEnd, i.isNegativeStrand, name)
        }.toJavaList
        val intervalList = new IntervalList(dict)
        intervalList.addall(intervals)
        intervalList.toList
    }
  }
}

/** Metric for [[AnnotateByIntervalList]].
  *
  * @param name             the name of the intervals
  * @param bases            the number of bases mapping to intervals with the given name.
  * @param unique_bases     the number of bases mapping exclusively to intervals with the given name,
  * @param span             the number of unique genomic loci/bases coverage by intervals with the given name.
  * @param frac_genome      the fraction of unique genomic loci/bases coverage by intervals with the given name.
  * @param frac_total_bases the fraction of mapped bases covered by intervals with the given name.
  * @param mean_coverage    the mean coverage of intervals with the given name.
  * @param templates        the number of templates mapping to the intervals with the given name
  * @param unique_templates the number of templates mapping exclusively to the intervals with the given name
  */
case class CoverageByIntervalMetrics(
  name: String,
  unique_bases: Long,
  bases: Long,
  span: Long,
  frac_genome: Double,
  frac_total_bases: Double,
  mean_coverage: Double,
  templates: Long,
  unique_templates: Long
) extends Metric

/**
  * Enum for determining how to associate a read that is enclosed in more than a single interval.
  * - Duplicate counting associates the template base with all enclosing intervals.
  * - MostSpecific counting associates the template bases with the smallest enclosing interval.
  * - LeastSpecific counting associates the template bases with the largest enclosing interval.
  */
sealed trait OverlapAssociationStrategy extends EnumEntry

object OverlapAssociationStrategy extends FgBioEnum[OverlapAssociationStrategy] {
  def values: IndexedSeq[OverlapAssociationStrategy] = findValues

  case object Duplicate extends OverlapAssociationStrategy
  case object MostSpecific extends OverlapAssociationStrategy
  case object LeastSpecific extends OverlapAssociationStrategy
}
