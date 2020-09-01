package com.fulcrumgenomics.rnaseq

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, IteratorToJavaCollectionsAdapter, javaIterableToIterator}
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{FilePath, PathToBam, PathToIntervals}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Metric, ProgressLogger}
import enumeratum.EnumEntry
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.util.{Interval, IntervalList, OverlapDetector}

/** Metric for [[AnnotateByIntervalList]].
  *
  * @param name the name of the intervals
  * @param bases the number of bases mapping to intervals with the given name.
  * @param span the number of unique genomic loci/bases coverage by intervals with the given name.
  * @param frac_genome the fraction of unique genomic loci/bases coverage by intervals with the given name.
  * @param frac_total_bases the fraction of mapped bases covered by intervals with the given name.
  * @param mean_coverage the mean coverage of intervals with the given name.
  */
case class CoverageByIntervalMetrics
( name: String,
  bases: Long,
  span: Long,
  frac_genome: Double,
  frac_total_bases: Double,
  mean_coverage: Double
) extends Metric

sealed trait OverlapAssociationStrategy extends EnumEntry
object OverlapAssociationStrategy extends FgBioEnum[OverlapAssociationStrategy] {
  def values: IndexedSeq[OverlapAssociationStrategy] = findValues
  case object AllOrNothing extends OverlapAssociationStrategy
  case object Factional extends OverlapAssociationStrategy
}

@clp(group=ClpGroups.Fastq, description=
  """
    |Get percentage of bases falling in each interval type.
  """)
class AnnotateByIntervalList
(
  @arg(flag='i', doc="Mapped BAM file.") val input: PathToBam,
  @arg(flag='l', doc="Path to intervals") val intervals: PathToIntervals,
  @arg(flag='o', doc="Path to the output metrics file") val output: FilePath,
  @arg(flag='R', doc="The regular expression to use to classify the name of the interval") val nameRegex: Option[String] = None,
  @arg(flag='u', doc="The name of the intervals not in the given set of intervals") val unmatched: String = "unmatched",
  @arg(flag='m', doc="The minimum mapping quality for a read to be included.") val minMappingQuality: Int = 20,
  @arg(flag='q', doc="The minimum base quality for a base to be included.") val minBaseQuality: Int = 0,
  @arg(flag='a', doc="The strategy to determine which interval overlapping templates should be associated with")
    val overlapAssociationStrategy: OverlapAssociationStrategy = OverlapAssociationStrategy.AllOrNothing
) extends FgBioTool with LazyLogging {

  Io.assertReadable(this.input)
  Io.assertReadable(this.intervals)
  Io.assertCanWriteFile(this.output)

  override def execute(): Unit = {
    val progress = ProgressLogger(logger, verb="Processed", noun="loci", unit=1000000)
    val reader = SamSource(this.input)
    val dict = reader.dict

    val (counter, overlapDetector) = {
      val _counter  = new SimpleCounter[String]()
      val intervals = buildIntervals(dict.asSam)
      intervals.foreach { interval => _counter.count(interval.getName, 0) }
      val _overlapDetector = OverlapDetector.create[Interval](intervals.iterator.toJavaList)
      (_counter, _overlapDetector)
    }

    //Templates are assigned to an interval together. Iterate through read one reads, but count bases from the mate.
    val iterator = reader.iterator
    val records = iterator.filterNot(record => record.duplicate
        || record.secondOfPair
        || record.secondary
        || record.supplementary
        || !record.pf
        || record.mapq < minMappingQuality)

    var numLociMultiMatched: Long = 0

    records.foreach { record: SamRecord =>
      val recordInterval    = new Interval(record.refName, record.start, record.end)
      val maybeMateInterval = record.mateEnd.map(end => new Interval(record.mateRefName, record.mateStart, end))
      val overlaps          = overlapDetector.getOverlaps(recordInterval).toSeq
      val maybeMateOverlaps = maybeMateInterval.map(overlapDetector.getOverlaps(_).toSeq)

      if (overlaps.isEmpty && maybeMateOverlaps.isEmpty) {
        counter.count(this.unmatched, record.length)
      } else {
        val names:Seq[String] = (overlaps.map(_.getName)
          ++ maybeMateOverlaps.map(_.map(_.getName)).getOrElse(Seq.empty[String])).distinct
        // If there is an overlap use the association strategy to determine which bases are associated with
        // which interval.
        if (names.length > 1) {
          //For intervals that only overlap one interval the bases are associated with that one.

          //All or nothing associates all bases of the read with one interval or the other.

          //Fractional associates only the bases overlapping the interval

          numLociMultiMatched += 1
        } else {
          // When there is only one interval matched we associate all intersected bases with the single interval
          names.foreach { name =>
            counter.count(name, overlaps.filter(_.getName == name).map(_.getIntersectionLength(recordInterval)).sum)
            maybeMateOverlaps.map(mateOverlaps =>
              counter.count(name, mateOverlaps.filter(_.getName == name).map(_.getIntersectionLength(maybeMateInterval.get)).sum))
          }
        }
      }
      progress.record(record.refName, record.start)
    }
    reader.close()

    logger.info(f"Found $numLociMultiMatched%,d loci with multiple classifications.")

    // Get the total # of bases covered by each set of intervals
    val lengthByIntervalName = new SimpleCounter[String]()
    overlapDetector.getAll.toSeq.groupBy(_.getName).foreach { case (name, intvs) =>
      val intervalList = new IntervalList(dict.asSam)
      intervalList.addall(intvs.iterator.toJavaList)
      val length = intervalList.uniqued().map(_.length()).sum
      lengthByIntervalName.count(name, length)
    }
    lengthByIntervalName.count(this.unmatched, dict.length - lengthByIntervalName.total)

    val totalBases = counter.total.toDouble
    val metrics = counter.map { case (name, bases) =>
      // The total # of genomic bases covered by the intervals with the given name
      val span = lengthByIntervalName.countOf(name)

      CoverageByIntervalMetrics(
        name             = name,
        bases            = bases,
        span             = span,
        frac_genome      = span / dict.length.toDouble,
        frac_total_bases = if (totalBases > 0) bases / totalBases else 0f,
        mean_coverage    = if (span > 0) bases / span.toDouble else 0f
      )
    }.toSeq.sortBy(_.name)
    Metric.write(this.output, metrics)
  }

  private def buildIntervals(dict: SAMSequenceDictionary): List[Interval] = {
    val inputIntervalList = IntervalList.fromPath(this.intervals)
    this.nameRegex.map(_.r) match {
      case None => inputIntervalList.uniqued(true).toList
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