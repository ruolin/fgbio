package com.fulcrumgenomics.rnaseq

import java.util

import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{FilePath, PathToBam, PathToIntervals}
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import htsjdk.samtools.util.{Interval, IntervalList, OverlapDetector, SamLocusIterator}
import com.fulcrumgenomics.FgBioDef.IteratorToJavaCollectionsAdapter
import com.fulcrumgenomics.FgBioDef.javaIterableToIterator
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.util.{Metric, ProgressLogger}
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.filter.{DuplicateReadFilter, FailsVendorReadQualityFilter, SamRecordFilter, SecondaryOrSupplementaryFilter}

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
  @arg(flag='q', doc="The minimum base quality for a base to be included.") val minBaseQuality: Int = 0
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
      val intervals = buildIntervals(dict)
      intervals.foreach { intv => _counter.count(intv.getName, 0) }
      val _overlapDetector = OverlapDetector.create[Interval](intervals.toIterator.toJavaList)
      (_counter, _overlapDetector)
    }
    val iterator = {
      val iter = new SamLocusIterator(reader.toSamReader)
      val filters = new util.ArrayList[SamRecordFilter]()
      filters.add(new SecondaryOrSupplementaryFilter)
      filters.add(new FailsVendorReadQualityFilter)
      filters.add(new DuplicateReadFilter)
      iter.setSamFilters(filters)
      iter.setMappingQualityScoreCutoff(this.minMappingQuality)
      iter.setQualityScoreCutoff(this.minBaseQuality)
      iter.setEmitUncoveredLoci(false)
      iter.setIncludeIndels(false)
      iter
    }

    var numLociMultiMatched: Long = 0

    iterator.foreach { locus: SamLocusIterator.LocusInfo =>
      val coverage          = locus.size
      val locusInterval     = new Interval(locus.getSequenceName, locus.getPosition, locus.getPosition)
      val overlaps          = overlapDetector.getOverlaps(locusInterval).toSeq
      if (overlaps.isEmpty) {
        counter.count(this.unmatched, coverage)
      } else {
        val names = overlaps.map(_.getName).distinct 
        if (names.length > 1) {
          logger.warning(
            s"Found more than one interval that matched locus '${locus.getSequenceName}:${locus.getPosition}' : ${names.mkString(", ")}"
          )
          numLociMultiMatched += 1
        }
        names.foreach { name => counter.count(name, coverage) }
      }
      progress.record(locus.getSequenceName, locus.getPosition)
    }
    reader.close()

    logger.info(f"Found $numLociMultiMatched%,d loci with multiple classifications.")

    // Get the total # of bases covered by each set of intervals
    val lengthByIntervalName = new SimpleCounter[String]()
    overlapDetector.getAll.toSeq.groupBy(_.getName).foreach { case (name, intvs) =>
      val intervalList = new IntervalList(dict)
      intervalList.addall(intvs.toIterator.toJavaList)
      val length = intervalList.uniqued().map(_.length()).sum
      lengthByIntervalName.count(name, length)
    }
    lengthByIntervalName.count(this.unmatched, dict.getReferenceLength - lengthByIntervalName.total)

    val totalBases = counter.total.toDouble
    val metrics = counter.map { case (name, bases) =>
      // The total # of genomic bases covered by the intervals with the given name
      val span = lengthByIntervalName.countOf(name)

      CoverageByIntervalMetrics(
        name             = name,
        bases            = bases,
        span             = span,
        frac_genome      = span / dict.getReferenceLength.toDouble,
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
        val intvs = inputIntervalList.map { i =>
          val name = i.getName match {
            case regex(n) => n
            case _        => throw new IllegalArgumentException(s"Could not match name '${i.getName}' with regex '${this.nameRegex}'")
          }
          new Interval(i.getContig, i.getStart, i.getEnd, i.isNegativeStrand, name)
        }.toJavaList
        val intervalList = new IntervalList(dict)
        intervalList.addall(intvs)
        intervalList.toList
    }
  }
}
