/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.vcf

import java.nio.file.Files
import java.util
import java.util.NoSuchElementException

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import dagr.commons.CommonsDef._
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt._
import htsjdk.samtools.util.StringUtil
import htsjdk.variant.variantcontext._
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.vcf._

import scala.collection.JavaConversions._
import scala.collection.JavaConverters._

@clp(
  description =
    """
      |Converts the output of HapCut to a VCF.
      |
      |The output of HAPCUT does not include all input variants, but simply those variants that are in phased blocks.
      |This tool takes the original VCF and the output of HapCut, and produces a VCF containing both the variants that
      |were phased and the variants that were not phased, such that all variants in the original file are in the output.
      |
      |The original VCF provided to HAPCUT is assumed to contain a single sample, as HAPCUT only supports a single
      |sample.
      |
      |By default, all phased genotypes are annotated with the "PS" (phase set) FORMAT tag, which by convention is the
      |position of the first variant in the phase set (see the VCF specification).  Furthermore, this tool formats the
      |alleles of a phased genotype using the '|' separator instead of the '/' separator, where the latter indicates the
      |genotype is unphased.  If the option to output phased variants in GATK's ReadBackedPhasing format is used, then
      |the first variant in a phase set will have '/' instead of '|' separating its alleles.  Also, all phased variants
      |will have "PASS" set in its FILTER column, while unphased variants will have "NotPhased" set in their FILTER
      |column.  Unlike GATK's ReadBackedPhasing, homozygous variants will always be unphased.
      |
      |More information about the purpose and operation of GATK's Read-backed phasing, including its output format, can
      |be found here:
      |  http://gatkforums.broadinstitute.org/gatk/discussion/45/purpose-and-operation-of-read-backed-phasing
      |
      |Additional FORMAT fields for phased variants are provided corresponding to per-genotype information produced by
      |HAPCUT:
      |  1. The "RC" tag gives the counts of calls supporting allele0 and allele1 respectively.
      |  2. The "LC" tag gives the change in likelihood if this SNP is made homozygous or removed.
      |  3. The "MCL" tag gives the maximum change in likelihood if this SNP is made homozygous or removed.
      |  4. The "RMEC" tag gives the reduction in MEC score if we remove this variant altogether.
      |
      |For more information about HAPCUT, see the source code or paper below.
      |  source code: https://github.com/vibansal/hapcut
      |  HapCut paper: An efficient and accurate algorithm for the haplotype assembly problem Bioinformatics. 2008 Aug
      |    15;24(16):i153-9.
      |
    """,
  group=ClpGroups.VcfOrBcf
)
class HapCutToVcf
( @arg(flag="v", doc="The original VCF provided to HAPCUT.") val vcf: PathToVcf,
  @arg(flag="i", doc="The output file from HAPCUT.") val input: FilePath,
  @arg(flag="o", doc="The output VCF with both phased and unphased variants.") val output: PathToVcf,
  @arg(flag="r", doc="Output phased variants in GATK's ReadBackedPhasing format.") val gatkPhasingFormat: Boolean = false
) extends FgBioTool with LazyLogging {
  Io.assertReadable(vcf)
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val vcfReader      = new VCFFileReader(vcf.toFile, false)
    val vcfWriter      = makeWriter(output, vcfReader)

    new HapCutAndVcfMergingIterator(input, vcfReader, gatkPhasingFormat).foreach(vcfWriter.add)

    vcfReader.safelyClose()
    vcfWriter.close()
  }

  /** Creates a VCF writer, adding extra header lines if the output is the phased VCF. */
  private def makeWriter(path: PathToVcf, vcfReader: VCFFileReader): VariantContextWriter = {
    val inputHeader = vcfReader.getFileHeader
    val builder = new VariantContextWriterBuilder()
      .setOutputFile(path.toFile)
      .setReferenceDictionary(inputHeader.getSequenceDictionary)
      .setOption(Options.INDEX_ON_THE_FLY)
      .modifyOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER, true)
    val writer: VariantContextWriter = builder.build

    // get the header lines in the input header that we wish to skip/replace with our own definitions
    val headerLinesToSkip = HapCutVcfHeaderLines.formatHeaderKeys.flatMap(key => Option(inputHeader.getFormatHeaderLine(key)))
    val headerLines: util.Set[VCFHeaderLine] = new util.HashSet[VCFHeaderLine](
      inputHeader.getMetaDataInSortedOrder.filterNot(headerLinesToSkip.contains).toSet
    )

    // add standard header lines
    VCFStandardHeaderLines.addStandardFormatLines(headerLines, false, Genotype.PRIMARY_KEYS)

    // add the new format header lines
    headerLines.addAll(HapCutVcfHeaderLines.formatHeaderLines)

    // add the new filter header line if we are to set a filter (not phased) on the unphased variants
    if (gatkPhasingFormat) headerLines.add(HapCutVcfHeaderLines.NotPhasedFilterHeaderLine)

    // create it
    val outHeader = new VCFHeader(headerLines, inputHeader.getSampleNamesInOrder) // create the header

    // write it and return
    writer.writeHeader(outHeader)
    writer
  }
}

object HapCutVcfHeaderLines {

  val ReadCountFormatTag = "RC"
  val ReadCountFormatDescription = "Counts of calls supporting allele0 and allele1 respectively"
  val ReadCountFormatHeaderLine = new VCFFormatHeaderLine(ReadCountFormatTag, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, ReadCountFormatDescription)

  val LikelihoodChangeFormatTag = "LC"
  val LikelihoodChangeFormatDescription = "Change in likelihood if this SNP is made homozygous or removed"
  val LikelihoodChangeFormatHeaderLine = new VCFFormatHeaderLine(LikelihoodChangeFormatTag, 3, VCFHeaderLineType.Float, LikelihoodChangeFormatDescription)

  val MaxLikelihoodChangeFormatTag = "MLC"
  val MaxLikelihoodChangeFormatDescription = "Maximum change in likelihood if this SNP is made homozygous or removed"
  val MaxLikelihoodChangeFormatHeaderLine = new VCFFormatHeaderLine(MaxLikelihoodChangeFormatTag, 1, VCFHeaderLineType.Float, MaxLikelihoodChangeFormatDescription)

  val MecReductionFormatTag = "RMEC"
  val MecReductionFormatDescription = "Reduction in MEC score if we remove this variant altogether"
  val MecReductionFormatHeaderLine = new VCFFormatHeaderLine(MecReductionFormatTag, 1, VCFHeaderLineType.Float, MecReductionFormatDescription)

  val PhaseSetFormatTag = "PS"
  val PhaseSetFormatDescription = "Phase set for the genotype (position of the first variant)"
  val PhaseSetFormatHeaderLine = new VCFFormatHeaderLine(PhaseSetFormatTag, 1, VCFHeaderLineType.Integer, PhaseSetFormatDescription)

  val NotPhasedFilterName        = "NotPhased"
  val NotPhasedFilterDescription = "The variant was not phased by HapCut."
  val NotPhasedFilterHeaderLine  = new VCFFilterHeaderLine(NotPhasedFilterName, NotPhasedFilterDescription)

  /** The VCF header FORMAT keys that will be added for HapCut specific-genotype information. */
  val formatHeaderKeys: Seq[String] = {
    Seq(ReadCountFormatTag, LikelihoodChangeFormatTag, MaxLikelihoodChangeFormatTag, MecReductionFormatTag, PhaseSetFormatTag)
  }

  /** The VCF header FORMAT lines that will be added for HapCut specific-genotype information. */
  val formatHeaderLines: Seq[VCFHeaderLine] = {
    Seq(ReadCountFormatHeaderLine, LikelihoodChangeFormatHeaderLine, MaxLikelihoodChangeFormatHeaderLine, MecReductionFormatHeaderLine, PhaseSetFormatHeaderLine)
  }
}

/** An iterator over all variants in the provided VCFFileReader.
  *
  * @param hapCutFile the output file with phased variants produced by HapCut.
  * @param vcfReader the VCF file with the original variants provided to HapCut.
  * @param gatkPhasingFormat true to output in GATK's ReadBackedPhasing format, false if to use the recommendations in the VCF spec.
  */
private class HapCutAndVcfMergingIterator(hapCutFile: FilePath, vcfReader: VCFFileReader, gatkPhasingFormat: Boolean = false)
  extends Iterator[VariantContext] {

  private val sourceIterator = vcfReader.toStream.zipWithIndex.iterator
  private val hapCutIterator = new HapCutReader(hapCutFile)
  private val sampleName = vcfReader.getFileHeader.getSampleNamesInOrder.head

  def hasNext(): Boolean = sourceIterator.hasNext

  /** Returns either the phased variant with added HapCut-specific genotype information (Left), or the original variant (Right). */
  def next(): VariantContext = {
    if (!hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false")
    val (sourceContext, sourceOffset) = sourceIterator.next()
    if (!hapCutIterator.hasNext()) formatSourceContext(sourceContext)
    else {
      hapCutIterator.next() match {
        case None => formatSourceContext(sourceContext)
        case Some(hapCutCall) =>
          if (hapCutCall.offset != sourceOffset+1) throw new IllegalStateException("BUG: calls are out of order")
          val hapCutContext = hapCutCall.toVariantContext(sampleName)
          replaceGenotypes(source=sourceContext, genotype=hapCutContext.getGenotype(0), isPhased = !gatkPhasingFormat || !hapCutCall.firstInBlock)
      }
    }
  }

  /** Returns a new variant context with the phase unset (if it was set) and appropriately formatted (filtered as
    * not phased if `gatkPhasingFormat` is true).
    */
  private def formatSourceContext(sourceContext: VariantContext): VariantContext = {
    val hasPhasedGenotype = sourceContext.getGenotypes.exists(_.isPhased)
    val hasPhasingSetId   = sourceContext.getGenotypes.exists(_.hasExtendedAttribute(HapCutVcfHeaderLines.PhaseSetFormatTag))

    if (!hasPhasedGenotype && !hasPhasingSetId && !gatkPhasingFormat) sourceContext
    else {
      val builder = new VariantContextBuilder(sourceContext)
      if (hasPhasedGenotype || hasPhasingSetId) {
        // unset the phase and remove the phasing set ID if the input has phase set
        builder.genotypes(sourceContext.getGenotypes().map { g =>
          val builder = new GenotypeBuilder(g).phased(false)
          val attrs = g.getExtendedAttributes.filterNot { case (tag, value) =>  tag == HapCutVcfHeaderLines.PhaseSetFormatTag }
          builder.noAttributes()
          builder.attributes(attrs)
          builder.make()
        })
      }
      if (gatkPhasingFormat) {
        // set the variant as filtered due to not being phased
        builder.filter(HapCutVcfHeaderLines.NotPhasedFilterName)
      }
      builder.make()
    }
  }

  /** Replaces the original genotype with HapCut's genotype and adds any HapCut-specific genotype information. */
  private def replaceGenotypes(source: VariantContext, genotype: Genotype, isPhased: Boolean): VariantContext = {
    val builder = new VariantContextBuilder(source)
    val sourceAlleles = source.getAlleles.toSeq
    val genotypeAlleles = genotype.getAlleles.map {
      allele => sourceAlleles.find(a => a.toString == allele.toString) match {
        case None =>
          throw new IllegalStateException(s"Could not find allele '$allele' in source alleles: " + sourceAlleles.map{_.toString}.mkString(", "))
        case Some(a) => a
      }
    }
    val sourceGenotype = source.getGenotype(0)
    val genotypeBuilder = new GenotypeBuilder(sourceGenotype).alleles(genotypeAlleles).phased(isPhased)
    genotypeBuilder.attributes(genotype.getExtendedAttributes)
    builder.genotypes(genotypeBuilder.make()).make()
  }
}

/** Haplotype block information produced by HapCut.
  *
  * @param phaseSetOption the phase set identifer, or None if not yet available.
  * @param offset the 1-based variant offset from the first variant given to HapCut.
  * @param len the total number of variants that this block spans, including unphased variants.
  * @param phased the total number of variants phased.
  * @param span the number of bases this block spans.
  * @param mecScore the minimum error correction (MEC) score.
  * @param fragments the number of fragments contributing to variants in this block.
  */
private case class HapCutBlockInfo (private val phaseSetOption: Option[Int] = None,
                                    offset: Int,
                                    len: Int,
                                    phased: Int,
                                    span: Int,
                                    mecScore: Float,
                                    fragments: Int) {
  /** The 1-based variant offset (from the first variant given to HapCut) of the last genotype in this block. */
  def offsetOfLastVariant: Int = this.offset + this.len - 1

  /** The phase set identifier for this block (typically the position of the first variant). Throws an exception if
    * the `phaseSetOption` was None when this was constructed. */
  def phaseSet: Int = phaseSetOption.get
}

/** Genotype-level information produced by HapCut */
private class HapCutCallGenotypeInfo private(val readCounts: List[Int], val likelihoods: List[Float], val delta: Float, val rMEC: Float) {
  import HapCutVcfHeaderLines._

  /** Adds the HapCut-specific genotype information to the given genotype builder. */
  def addTo(builder: GenotypeBuilder): GenotypeBuilder = {
    builder
      .attribute(ReadCountFormatTag, this.readCounts.asJava)
      .attribute(LikelihoodChangeFormatTag, this.likelihoods.asJava)
      .attribute(MaxLikelihoodChangeFormatTag, this.delta)
      .attribute(MecReductionFormatTag, this.rMEC)
  }
}

/** Values and methods for HapCut-specific genotype information. */
private object HapCutCallGenotypeInfo {
  // A developer's notes
  //  R0,R1:L00,L01,L11:delta:rMEC
  //  "%d,%d:%0.1f,%0.1f,%0.1f:%0.1f:%0.1f"
  //  R0,R1: counts of calls supporting allele0 and allele1 respectively
  //  L00,L01,L11: change in likelihood if this SNP is made homozygous or removed
  //  delta: maximum change in likelihood if this SNP is made homozygous or removed
  //  rMEC: reduction in MEC score if we remove this variant altogether
  val HapCutCallInfoPattern = """(\d+),(\d+):(\d+\.\d+),(\d+\.\d+),(\d+\.\d+):(\d+\.\d+):(\d+\.\d+)""".r

  /** Parses the genotype information produced by HapCut. */
  def apply(info: String): HapCutCallGenotypeInfo = {
    val tokens = info.split(":")
    new HapCutCallGenotypeInfo(
      readCounts  = tokens(0).split(",").map(_.toInt).toList,
      likelihoods = tokens(1).split(",").map(_.toFloat).toList,
      delta       = tokens(2).toFloat,
      rMEC        = tokens(3).toFloat
    )
  }
}

/** Information about a variant phased by HapCut.
  *
  * @param offset the 1-based offset of the variant relative to the variants given to HapCut.
  * @param hap1Allele the first allele.
  * @param hap2Allele the second allele.
  * @param contig the contig.
  * @param pos the position.
  * @param ref the reference allele.
  * @param alts the alternate alleles.
  * @param genotype the genotype.
  * @param info the HapCut-specific genotype information.
  * @param phaseSet the phase set identifier for the phase block to which this call belongs
  */
private case class HapCutCall private(offset: Int,
                                      hap1Allele: Int, hap2Allele: Int,
                                      contig: String, pos: Int,
                                      ref: String, alts: String,
                                      genotype: String, info: HapCutCallGenotypeInfo,
                                      phaseSet: Int
                                     ) {
  import HapCutVcfHeaderLines._

  /** true if this variant is the first phased variant in a phased block, false otherwise. */
  def firstInBlock: Boolean = phaseSet == pos

  /** Converts the HapCut variant representation to a VariantContext. */
  def toVariantContext(sampleName: String): VariantContext = {
    // Parse the alleles
    val refAllele = Allele.create(this.ref, true)
    val altAlleles = this.alts.split(",").map { alt =>
      Allele.create(alt, false)
    }.toSeq
    val alleles = refAllele +: altAlleles
    val allelesCollection = alleles.asJavaCollection

    // Assumes the genotype is before the first ":" token
    val genotypeAlleles = this.genotype.split(":").head.split("/").map(_.toInt).map(i => alleles(i))
    val genotype = this.info.addTo(new GenotypeBuilder(sampleName, genotypeAlleles.toList))
      .attribute(PhaseSetFormatTag, phaseSet)
      .make()

    // build the context, and make sure to recompute the end
    new VariantContextBuilder(
      s"${this.offset}",
      this.contig,
      this.pos.toLong,
      this.pos.toLong,
      allelesCollection
    )
      .computeEndFromAlleles(seqAsJavaList(alleles), this.pos)
      .genotypes(genotype).make()
  }
}

private object HapCutCall {

  /** Parse a variant line.
    *
    * @param callLine the line to parse.
    * @param firstVariantPosition the position of the first variant in the block to which this variant belongs, or None
    *                             if this variant is the first variant in the block.
    **/
  def apply(callLine: String, firstVariantPosition: Option[Int]): HapCutCall = {
    val tokens = new Array[String](9)
    val numTokens = StringUtil.split(callLine, tokens, '\t')
    if (numTokens != tokens.length) {
      throw new IllegalStateException(s"Did not parse 9 fields (parsed $numTokens) in the HapCut call: " + callLine.trim)
    }
    val offset = tokens(0).toInt
    val pos    = tokens(4).toInt
    new HapCutCall(
      offset       = offset,
      hap1Allele   = tokens(1).toInt,
      hap2Allele   = tokens(2).toInt,
      contig       = tokens(3),
      pos          = pos,
      ref          = tokens(5),
      alts         = tokens(6),
      genotype     = tokens(7),
      info         = HapCutCallGenotypeInfo(tokens(8)),
      phaseSet     = firstVariantPosition.getOrElse(pos)
    )
  }
}

object HapCutReader {
  val BlockLinePattern = """^BLOCK: offset: (\d+) len: (\d+) phased: (\d+) SPAN: (\d+) MECscore (\d+.\d+) fragments (\d+)$""".r
}

/** Reads phased blocks and variants from a HapCut output file.  For each variant call given to HapCut, returns either
  * some `HapCutCall` if that variant was phased, or None if it was not phased.  HapCut does not include unphased
  * variants in its output, but does provide the 1-based offset of every phased variant relative to the variant calls
  * given to HapCut.  Therefore, we are able to determine if None should be returned (unphased, not in the HapCut file)
  * or if a `HapCutCall` should be returned (phased, in the HapCut file).  Nonetheless, variants that were not phased
  * after the last HapCut phased variant will not be returned (hasNext will be false), as HapCut doesn't provide the
  * total # of variants that were given to HapCut, so these trailing unphased variant calls must be treated specially
  * for now.
  */
private[vcf] class HapCutReader(input: FilePath) extends Iterator[Option[HapCutCall]] {
  private val lineStream = Files.lines(input)
  private val lineIterator = lineStream.iterator().buffered
  private var block: Option[HapCutBlockInfo] = None
  private var calls: Iterator[Option[HapCutCall]] = Iterator.empty

  // initialize
  this.advance()

  def hasNext(): Boolean = {
    this.block match {
      case Some(b) if calls.hasNext => true
      case _ =>
        if (this.lineIterator.hasNext) this.advance()
        else false
    }
  }

  def next(): Option[HapCutCall] = {
    if (!hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false")
    calls.next()
  }

  def close(): Unit = {
    this.lineStream.close()
  }

  private def nextBlock(): Unit = {
    // Read the block information
    val blockLine    = this.lineIterator.next().trim
    val currentBlock = blockLine match {
      case HapCutReader.BlockLinePattern(offset, len, phased, span, mecScore, fragments) =>
        new HapCutBlockInfo(None, offset.toInt, len.toInt, phased.toInt, span.toInt, mecScore.toFloat, fragments.toInt)
      case _ => throw new IllegalStateException("Could not parse block line: " + blockLine)
    }

    // Read in the phased variants, storing a map of their offset (1-based) in the original file to the HapCutCall.
    var firstVariantInBlock: Option[Int] = None
    val callMap = (for (i <- 1 to currentBlock.phased) yield {
      if (!lineIterator.hasNext) throw new IllegalStateException("Reached the end of file when searching for a call")
      val call = HapCutCall(lineIterator.next(), firstVariantInBlock)
      if (firstVariantInBlock.isEmpty) firstVariantInBlock = Some(call.pos)
      (call.offset, call)
    }).toMap

    // add in any empty calls
    val nextVariantOffset = this.block.map(_.offsetOfLastVariant).getOrElse(0) + 1
    this.calls = Stream.range(nextVariantOffset, currentBlock.len+currentBlock.offset).map(callMap.get).toIterator.buffered

    // update the current block
    this.block = Some(currentBlock.copy(phaseSetOption=firstVariantInBlock))
  }

  /** Reads in the next block of data */
  private def advance(): Boolean = {
    // find the next block
    while (this.lineIterator.hasNext && !this.lineIterator.head.startsWith("BLOCK:")) {
      this.lineIterator.next()
    }
    if (this.lineIterator.hasNext) { nextBlock(); true }
    else false
  }
}
