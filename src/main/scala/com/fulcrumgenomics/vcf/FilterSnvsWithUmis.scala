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
import com.fulcrumgenomics.vcf.FilterSnvsWithUmis.MidPileup
import dagr.commons.CommonsDef.{PathToBam => _, PathToVcf => _, yieldAndThen => _}
import dagr.sopt.{arg, clp}
import htsjdk.samtools.util.SamLocusIterator.LocusInfo
import htsjdk.samtools.util.{Interval, IntervalList, SamLocusIterator, SequenceUtil}
import htsjdk.samtools.{SamReader, SamReaderFactory}
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
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
  @arg(flag="e", doc="A floor on the computed error rate.") val errorRate: Double = 1/1e6,
  @arg(flag="c", doc="Width of condfidence interval for error rate computation.") val confidence: Double = 0.95
// TODO: base quality and mapq cutoffs?
// TODO: false positive rate?
) extends FgBioTool {



  override def execute(): Unit = {
    val in  = new VCFFileReader(input.toFile, false)
    // TODO add header lines for a) new filter, b) a new INFO field to store the computed error rate
    val out = new VariantContextWriterBuilder().setOutputFile(output.toFile).setOption(Options.INDEX_ON_THE_FLY).build()
    val sam = SamReaderFactory.make().open(bam)

    in.foreach { variant =>
      if (variant.isSNP) {
        val gt  = variant.getGenotype(0)
        if (gt.isHet) {
          val a1 = gt.getAllele(0).getBases()(0)
          val a2 = gt.getAllele(1).getBases()(0)
          val ad  = gt.getAD
          val af  = min(ad(0), ad(1)) / (ad(0) + ad(1)).toDouble
          val pile = pileup(sam, variant.getContig, variant.getStart)
          val midPiles = reduce(pile)
          val (rate, depth) = calculateErrorRateAndDepth(midPiles, a1, a2)

//          if (looksLikeFalsePositive(depth=depth, af = ???, errorRate=rate, 1/1e6))
        }
      }
      else {
        out.add(variant)
      }
    }
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
    iterator.setMappingQualityScoreCutoff(0)
    iterator.setQualityScoreCutoff(0)
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
      val base = rec.getReadBase
      if (!names.contains(name)) {
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
  def calculateErrorRateAndDepth(piles: Seq[MidPileup], allele1: Byte, allele2: Byte): (Double, Int) = {
    val count1 = piles.map(p => p(allele1)).sum
    val count2 = piles.map(p => p(allele2)).sum

    val (major, minor) = if (count1 > count2) (allele1, allele2) else (allele2, allele1)
    var (total, errors) = (0,0)

    piles.foreach { pile =>
      val majorCount = pile(major)
      val minorCount = pile(minor)
      total += majorCount + minorCount
      if (majorCount > 0 && minorCount > 0) errors += minorCount
    }

    (Math.max(this.errorRate, errors / total.toDouble), total)
  }

  def looksLikeFalsePositive(depth: Int, af: Double, errorRate: Double, pvalue: Double): Boolean = {
    val dist   = new BinomialDistribution(depth, errorRate)
    val altObs = Math.floor(depth * af).toInt
    val pError = 1 - dist.cumulativeProbability(altObs - 1)
    pError < pvalue
  }
}

