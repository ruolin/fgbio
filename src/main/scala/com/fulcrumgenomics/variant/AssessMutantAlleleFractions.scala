package com.fulcrumgenomics.variant

import java.nio.file.Path

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.Io
import dagr.commons.CommonsDef.{PathToBam, PathToIntervals, PathToVcf}
import dagr.sopt.{arg, clp}
import htsjdk.samtools.filter.DuplicateReadFilter
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset
import htsjdk.samtools.util._
import htsjdk.samtools.{SAMFileHeader, SAMSequenceDictionary, SamReader, SamReaderFactory}
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.JavaConversions._
import scala.collection.mutable
import scala.collection.mutable.ListBuffer
import scala.util.Random

/**
  * Tool to check the mutant allele frequencies at a set of SNP sites described by a VCF. Each site
  * is downsampled (once) to the target coverage if it is above the specified coverage, and then the
  * reads are examined and counts of the ref and alt alleles produced.
  *
  * Output is a simple tab-separated text table with one row per SNP.
  */
@clp(
  description =
    """
      |Stuff goes here
    """,
  group=ClpGroups.VariantAssessment
)
class AssessMutantAlleleFractions
(
  @arg(flag="i", doc="The input SAM or BAM file.") val input: PathToBam,
  @arg(flag="v", doc="A VCF containing the set of sites to extract mutant allele fractions at.") val vcf: PathToVcf,
  @arg(flag="l", doc="A set of regions to restrict analysis to.") val intervals: Option[PathToIntervals] = None,
  @arg(flag="o", doc="The output file to write to.") val output: Path,
  @arg(flag="s", doc="Use the AF field sample's genotype for expected allele frequency.") val sample: Option[String] = None,
  @arg(flag="d", doc="The target coverage to downsample to.") val downsampleTo: Int = 1000,
  @arg(flag="m", doc="Exclude reads with mapping quality below this value.") val minMappingQuality: Int = 10,
  @arg(flag="q", doc="Exclude bases below this quality value.") val minQuality: Int = 10,
  @arg(doc="Allow bases to be counted from overlapping reads from the same insert.") val allowOverlappingReads: Boolean = false,
  @arg(doc="If true include duplicate reads, otherwise exclude duplicate reads.") val includeDuplicates: Boolean = false
) extends FgBioTool {

  override def execute(): Unit = {
    // Input checking
    Io.assertReadable(Seq(input, vcf))
    intervals.foreach(Io.assertReadable)
    Io.assertCanWriteFile(output)
    assert(downsampleTo > 0, "Downsampling target must be > 0.")

    // Open up all the files and objects needed to iterate over the data
    val in          = SamReaderFactory.make.open(input.toFile)
    val sites       = buildIntervalList(vcf, in.getFileHeader.getSequenceDictionary)
    val vcfIterator = new VCFFileReader(vcf.toFile, false).iterator.filter(_.isSNP)
    val bamSample   = in.getFileHeader.getReadGroups.iterator.next.getSample
    val fmt         = new FormatUtil
    val out         = Io.toWriter(output)


    out.append("chrom\tposition\tid\tref_allele\talt_allele\tsample\texpected_maf\tdepth\tref_count\talt_count\tother_count\tmaf\n")
    val samLocusIterator = buildLocusIterator(in, sites).iterator()

    for (info <- samLocusIterator) {
      val ctx: VariantContext = vcfIterator.next()
      if (info.getSequenceName != ctx.getContig || info.getPosition != ctx.getStart) {
        throw new IllegalStateException("BAM and VCF out of sync. Bam @ " + info.toString + ", VCF @ " + ctx.getContig + ":" + ctx.getStart)
      }
      val refAllele: Byte = ctx.getReference.getBases.apply(0)
      val altAllele: Byte = ctx.getAlternateAllele(0).getBases.apply(0)
      val expectedAf: String = sample match {
        case Some(s) =>
          Option(ctx.getGenotype(s).getAnyAttribute("AF")).map(_.asInstanceOf[String]).getOrElse("")
        case None if ctx.getNSamples == 1 =>
          Option(ctx.getGenotype(ctx.getSampleNames.head).getAnyAttribute("AF")).map(_.asInstanceOf[String]).getOrElse("")
        case None => ""
      }

      var (refCount, altCount, otherCount) = (0,0,0)
      downsample(filterForOverlaps(info.getRecordAndPositions), downsampleTo).foreach(_.getReadBase match {
        case `refAllele` => refCount += 1
        case `altAllele` => altCount += 1
        case _           => otherCount += 1
      })

     val maf: Double = altCount / (altCount + refCount + otherCount).toDouble
      out.append(info.getSequenceName).append('\t')
      out.append(fmt.format(info.getPosition)).append('\t')
      out.append(ctx.getID).append('\t')
      out.append(refAllele.toChar).append('\t')
      out.append(altAllele.toChar).append('\t')
      out.append(bamSample).append('\t')
      out.append(expectedAf).append('\t')
      out.append(fmt.format(refCount + altCount + otherCount)).append('\t')
      out.append(fmt.format(refCount)).append('\t')
      out.append(fmt.format(altCount)).append('\t')
      out.append(fmt.format(otherCount)).append('\t')
      out.append(fmt.format(maf))
      out.append('\n')
    }
    out.close
  }

  /**
    * Ensures that the returned list only contains one record per read-name, effectively
    * eliminating any double-counting where reads from the same template overlap the same position.
    */
  private def filterForOverlaps(records: TraversableOnce[RecordAndOffset]): Seq[RecordAndOffset] = {
    if (this.allowOverlappingReads) {
      records.toSeq
    }
    else {
      val shuffled = Random.shuffle(records)
      val retval = ListBuffer[RecordAndOffset]()
      val names = mutable.Set[String]()
      shuffled.foreach(rec => if (names.add(rec.getRecord.getReadName)) retval.add(rec))
      return retval.toList
    }
  }

  /** Builds a SamLocusIterator with appropriate filtering. */
  private def buildLocusIterator(in: SamReader, sites: IntervalList): SamLocusIterator = {
    val iterator: SamLocusIterator = new SamLocusIterator(in, sites)
    iterator.setSamFilters( if(includeDuplicates) Nil else List(new DuplicateReadFilter))
    iterator.setEmitUncoveredLoci(true)
    iterator.setMappingQualityScoreCutoff(minMappingQuality)
    iterator.setQualityScoreCutoff(minQuality)
    iterator
  }

  /** Builds an interval list out of the entries in a VCF. */
  private def buildIntervalList(vcf: PathToVcf, dict: SAMSequenceDictionary): IntervalList = {
    val in: VCFFileReader = new VCFFileReader(vcf.toFile, false)
    val header: SAMFileHeader = new SAMFileHeader
    header.setSequenceDictionary(dict)
    val ilist: IntervalList = new IntervalList(header)

    in.iterator().filter(_.isSNP)
      .map(vc => new Interval(vc.getContig, vc.getStart, vc.getEnd, false, vc.getID))
      .foreach(ilist.add)

    ilist.sorted
  }

  /**
    * Downsamples a list of RecordAndOffset to a given coverage. If the list is already below the requested
    * size, it is returned as is.
    */
  private[variant] def downsample[T](in: Seq[T], target: Int): Seq[T] = {
    if (in.size <= target) return in
    Random.shuffle(in).take(target)
  }
}

