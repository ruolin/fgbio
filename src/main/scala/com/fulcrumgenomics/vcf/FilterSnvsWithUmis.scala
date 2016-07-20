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

package com.fulcrumgenomics.vcf

import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.ProgressLogger
import com.fulcrumgenomics.vcf.FilterSnvsWithUmis.MidPileup
import dagr.commons.CommonsDef.{PathToBam => _, PathToVcf => _, yieldAndThen => _}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.util.SamLocusIterator.LocusInfo
import htsjdk.samtools.util.{Interval, IntervalList, SamLocusIterator, SequenceUtil}
import htsjdk.samtools.{SamReader, SamReaderFactory}
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.variantcontext.{VariantContext, VariantContextBuilder}
import htsjdk.variant.vcf._
import org.apache.commons.math3.distribution.BinomialDistribution

import scala.collection.JavaConversions.iterableAsScalaIterable
import scala.collection.mutable
import scala.math.min

object FilterSnvsWithUmis {
  case class MidPileup(molecularId: String, as: Int, cs: Int, gs: Int, ts: Int) {
    def depth: Int = as + cs + gs + ts

    def apply(base: Byte): Int = SequenceUtil.upperCase(base).toChar match {
      case 'A' => as
      case 'C' => cs
      case 'G' => gs
      case 'T' => ts
      case _   => 0
    }
  }
}

@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Filters the SNVs in a VCF by using reads grouped by Molecular ID in a BAM file to compute
    |a local error rate and then filter variants by that error rate.
    |
    |INDELS and complex events are passed through
    |USES AD FIELD
    |FILTERING OF READS/BASES
  """
)
class FilterSnvsWithUmis
( @arg(flag="i", doc="A VCF containing variant calls for a single sample.") val input: PathToVcf,
  @arg(flag="o", doc="Output VCF to write.") val output: PathToVcf,
  @arg(flag="b", doc="A BAM file containing reads tagged with molecular IDs.") val bam : PathToBam,
  @arg(flag="t", doc="The SAM/BAM tag containing the unique molecule ID after grouping.") val tag: String = "MI",
  @arg(flag="e", doc="A floor on the computed error rate.")                       val errorRate: Double = 1/1e6,
  @arg(flag="f", doc="Acceptable false positive rate, e.g. 1 per mb = 1/1e6.")    val falsePositiveRate: Double = 1/1e6,
  @arg(flag="q", doc="Base quality cutoff for bases used to compute error rate.") val quality: PhredScore = 5.toByte,
  @arg(flag="x", doc="Remove filters from input VCF prior to filtering.")         val resetFilters: Boolean = false
) extends FgBioTool with LazyLogging {

  private def relevant(v: VariantContext): Boolean = v.isSNP

  override def execute(): Unit = {
    val in  = new VCFFileReader(input.toFile, false)
    val out = makeWriter(output, in)
    val sam = SamReaderFactory.make().open(bam)
    val progress = new ProgressLogger(logger=logger, noun="variants", unit=50)

    for (variant <- in; gt <- Option(variant.getGenotype(0))) {
      if (variant.isSNP && gt.isHet) {
        val a1 = gt.getAllele(0).getBases()(0)
        val a2 = gt.getAllele(1).getBases()(0)
        val ad  = gt.getAD
        val depth = ad(0) + ad(1)
        val af  = min(ad(0), ad(1)) / depth.toDouble
        val pile = pileup(sam, variant.getContig, variant.getStart)
        val midPiles = reduce(pile)
        val rate = calculateErrorRate(midPiles, a1, a2)

        val builder = new VariantContextBuilder(variant)
        if (this.resetFilters) builder.unfiltered()
        builder.attribute("ER", rate)

        if (looksLikeFalsePositive(depth=depth, af=af, errorRate=rate, falsePositiveRate)) builder.filter("LikelyError")
        else if (this.resetFilters) builder.passFilters()

        out.add(builder.make())
      }
      else {
        out.add(variant)
      }

      progress.record(variant.getContig, variant.getStart)
    }

    out.close()
    in.safelyClose()
  }

  /** Builds the output VCF writer. */
  private def makeWriter(path: PathToVcf, vcfReader: VCFFileReader): VariantContextWriter = {
    val inputHeader = vcfReader.getFileHeader
    val builder = new VariantContextWriterBuilder()
      .setOutputFile(path.toFile)
      .setReferenceDictionary(inputHeader.getSequenceDictionary)
      .setOption(Options.INDEX_ON_THE_FLY)
      .modifyOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER, true)
    val writer: VariantContextWriter = builder.build

    val lines = new util.HashSet[VCFHeaderLine](inputHeader.getMetaDataInSortedOrder)
    lines.add(new VCFInfoHeaderLine("ER", 1, VCFHeaderLineType.Float, "The computed error rate at the site."))
    lines.add(new VCFFilterHeaderLine("LikelyError", "The variant is likely the result of the local error rate."))

    val outHeader = new VCFHeader(lines, inputHeader.getSampleNamesInOrder)
    writer.writeHeader(outHeader)
    writer
  }

  /** Constructs a pileup over a single genomic base. */
  def pileup(in: SamReader, chrom: String, pos: Int): LocusInfo = {
    val interval = new Interval(chrom, pos, pos)
    val ilist    = new IntervalList(in.getFileHeader)
    ilist.add(interval)
    val iterator = new SamLocusIterator(in, ilist)
    iterator.setEmitUncoveredLoci(true)
    iterator.setMappingQualityScoreCutoff(1)
    iterator.setQualityScoreCutoff(this.quality)
    iterator.setSamFilters(new util.ArrayList())
    iterator.setMaxReadsToAccumulatePerLocus(Int.MaxValue)

    yieldAndThen(iterator.iterator().next()) { iterator.close() }
  }

  /** Reduces a LocusInfo into a small object with an ID and a count of each base observation. */
  def reduce(info: LocusInfo): Seq[MidPileup] = info.getRecordAndPositions.groupBy(x => x.getRecord.getStringAttribute("MI")).map { case (mi, recs) =>
    val names = mutable.Set[String]()
    var (a,c,g,t) = (0,0,0,0)
    recs.foreach { rec =>
      val name = rec.getRecord.getReadName
      if (!names.contains(name)) {
        val base = rec.getReadBase
        names += name
        SequenceUtil.upperCase(base).toChar match {
          case 'A' => a += 1
          case 'C' => c += 1
          case 'G' => g += 1
          case 'T' => t += 1
          case _   => Unit
        }
      }
    }

    MidPileup(molecularId=mi, as=a, cs=c, gs=g, ts=t)
  }.toSeq

  /** Calculates the error rate between two bases given the molecular pileups. */
  def calculateErrorRate(piles: Seq[MidPileup], allele1: Byte, allele2: Byte): Double = {
    val count1 = piles.map(p => p(allele1) / p.depth.toDouble).sum
    val count2 = piles.map(p => p(allele2) / p.depth.toDouble).sum

    val (major, minor) = if (count1 > count2) (allele1, allele2) else (allele2, allele1)
    var (total, errors) = (0,0)

    piles.foreach { pile =>
      val majorCount = pile(major)
      val minorCount = pile(minor)
      total += majorCount + minorCount
      if (majorCount > 0 && minorCount > 0) errors += minorCount
    }

    Math.max(this.errorRate, errors / total.toDouble)
  }

  def looksLikeFalsePositive(depth: Int, af: Double, errorRate: Double, pvalue: Double): Boolean = {
    val dist   = new BinomialDistribution(null, depth, errorRate)
    val altObs = Math.floor(depth * af).toInt
    val pError = 1 - dist.cumulativeProbability(altObs - 1)
    pError >= pvalue
  }
}

