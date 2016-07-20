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
import java.util.Collections

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.vcf.MakeTwoSampleMixtureVcf._
import dagr.sopt
import dagr.sopt.{arg, clp}
import htsjdk.samtools.util.IntervalList
import htsjdk.variant.variantcontext._
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.vcf._

import scala.collection.JavaConversions.asJavaCollection
import scala.collection.mutable.ListBuffer

object MakeTwoSampleMixtureVcf {
  val AlleleFractionField = "AF"
  val GermlineFilter      = "alt_allele_in_normal"
  val UnknownGtFilter     = "unknown_gt"
}

/**
  * Creates a VCF by mixing two germline samples at a given proportion.
  *
  * @author Tim Fennell
  */
@sopt.clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Creates a VCF by in-silico mixing data from two samples. Writes a new VCF for a 'tumor' and
    |optionally a 'normal' where the tumor is a mixture of two germline samples, and the normal
    |is one of the same germline samples. Only outputs loci that are variant in the 'tumor'.
    |Writes the 'AF' attribute in the format field with the expected allele fraction in the tumor."
  """
)
class MakeTwoSampleMixtureVcf
( @arg(flag="i", doc="Input VCF file.")  val input: PathToVcf,
  @arg(flag="o", doc="Output VCF file.") val output: PathToVcf,
  @arg(flag="t", doc="Name of the 'tumor' sample in the input VCF.")  val tumor: String,
  @arg(flag="n", doc="Name of the 'normal' sample in the input VCF.") val normal: String,
  @arg(flag="f", doc="What fraction of the mixture comes from the 'tumor' sample.") val tumorFraction: Double = 0.5,
  @arg(flag="T", doc="Tumor only mode - only output tumor genotypes and don't filter sites.") val tumorOnly: Boolean = false,
  @arg(flag="N", doc="Treat no-calls for either sample as hom-ref genotypes.") val noCallIsHomRef: Boolean = true,
  @arg(flag="l", doc="Optional set of intervals to restrict to.") val intervals: Option[PathToIntervals] = None
) extends FgBioTool {
  if (tumorFraction < 0 || tumorFraction > 1) fail("Tumor fraction must be between 0 and 1.")
  if (tumor == normal) fail("Tumor and Normal samples must be different samples.")

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  intervals.foreach(Io.assertCanWriteFile(_))

  override def execute(): Unit = {
    val in = new VCFFileReader(input.toFile, true)
    if (!in.getFileHeader.getSampleNamesInOrder.contains(tumor))  fail(s"VCF does not contain tumor sample ${tumor}")
    if (!in.getFileHeader.getSampleNamesInOrder.contains(normal)) fail(s"VCF does not contain normal sample ${normal}")

    val out = if (tumorOnly) makeWriter(in.getFileHeader, "tumor") else makeWriter(in.getFileHeader, "tumor", "normal")
    iterator(in, this.intervals).filterNot(_.getAlternateAlleles.size() > 1).foreach { ctx =>
      val oldTumorGt  = ctx.getGenotype(tumor)
      val oldNormalGt = ctx.getGenotype(normal)
      val oldGts = Seq(oldTumorGt, oldNormalGt)

      val refAllele = ctx.getReference
      val altAllele = ctx.getAlternateAllele(0)

      // Build the Tumor genotype
      val newTumorGt = {
        val builder = new GenotypeBuilder("tumor")
        val alleles = new util.ArrayList[Allele]()
        if (noCallIsHomRef || oldGts.forall(_.isCalled)) {
          if (oldGts.map(_.countAllele(refAllele)).sum + oldGts.count(_.isNoCall) > 0) alleles.add(refAllele)
          if (oldGts.map(_.countAllele(altAllele)).sum > 0) alleles.add(altAllele)
          if (alleles.size() == 1) alleles.add(alleles.get(0))
          builder.alleles(alleles)

          // We count alts only, because ref and no-call can be treated as ref here
          val af = (this.tumorFraction * oldTumorGt.countAllele(altAllele) / 2) + ((1-this.tumorFraction) * oldNormalGt.countAllele(altAllele) / 2)
          builder.attribute(MakeTwoSampleMixtureVcf.AlleleFractionField, af)
        }
        else {
          builder.alleles(util.Arrays.asList(Allele.NO_CALL, Allele.NO_CALL))
        }

        builder.make()
      }

      // Build the normal genotype
      val newNormalGt = tumorOnly match {
        case true  => None
        case false =>
          val builder = new GenotypeBuilder("normal")
          builder.alleles(oldNormalGt.getAlleles)
          Some(builder.make())
      }

      val gts     = Seq(Some(newTumorGt), newNormalGt).flatten
      val builder = new VariantContextBuilder(ctx.getSource, ctx.getContig, ctx.getStart(), ctx.getEnd(), ctx.getAlleles)
      builder.genotypes(gts)

      val filters = ListBuffer[String]()
      if (!tumorOnly && !oldNormalGt.isHomRef) filters += GermlineFilter
      if (newTumorGt.isNoCall) filters += UnknownGtFilter
      if (filters.nonEmpty) builder.filters(filters:_*) else builder.passFilters()
      out.add(builder.make());
    }

    out.close()
    in.safelyClose()
  }

  /** Builds a header that can be used to write the mixture VCF. */
  def makeWriter(in: VCFHeader, samples: String*): VariantContextWriter = {
    import scala.collection.JavaConversions.{iterableAsScalaIterable}
    val header = new VCFHeader(Collections.emptySet[VCFHeaderLine](), util.Arrays.asList(samples:_*))

    in.getFilterLines.foreach(header.addMetaDataLine)
    header.addMetaDataLine(new VCFFilterHeaderLine(GermlineFilter,  "Evidence seen in the normal sample"))
    header.addMetaDataLine(new VCFFilterHeaderLine(UnknownGtFilter, "GT and AF not known due to one or more input samples being no-called."))

    in.getFormatHeaderLines.foreach(header.addMetaDataLine)
    header.addMetaDataLine(new VCFFormatHeaderLine("AF", 1, VCFHeaderLineType.Float, "Alt allele fraction in the tumor"))

    in.getInfoHeaderLines.foreach(header.addMetaDataLine)
    header.setSequenceDictionary(in.getSequenceDictionary)

    val writer = new VariantContextWriterBuilder()
      .setOutputFile(output.toFile)
      .setOption(Options.INDEX_ON_THE_FLY)
      .setReferenceDictionary(in.getSequenceDictionary)
      .build()

    writer.writeHeader(header)
    writer
  }

  /** Generates an iterator over the whole file, or over the intervals if provided. */
  def iterator(in: VCFFileReader, intervals: Option[PathToIntervals]): Iterator[VariantContext] = {
    import scala.collection.JavaConversions.asScalaIterator
    intervals match {
      case None       => in.iterator()
      case Some(path) => MultiIntervalVcfIterator(in, IntervalList.fromFile(path.toFile))
    }
  }
}
