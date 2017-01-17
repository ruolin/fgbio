/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.personal.tfenne

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, Metric}
import com.fulcrumgenomics.vcf.JointVariantContextIterator
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.JavaConversions.iterableAsScalaIterable

case class VariantInfo
( chrom: String,
  chr: String,
  pos: Int,
  ref: String,
  alt: String,
  mutation: String,
  variant_type: String,
  qual: Option[Double],
  filter: Option[String],
  conc_state: String,
  conc_desc: String,
  expected_af: Option[Double],
  genotype: String,
  ref_obs: Option[Int],
  alt_obs: Option[Int],
  af: Option[Double]
) extends Metric

@clp(group=ClpGroups.Personal, description=
  """
    |Creates a more human-readable tabular text file from genotype concordance data.  Takes as input
    |both the VCF produced by running Picard's GenotypeConcordance tool with OUTPUT_VCF=true and the
    |original call VCF.  Produces a text file with one line per variant with key information on each
    |variant designed to make filtering, grouping and plotting easier in R or Excel.
    |
    |NOTE: This tool exists in the 'personal' group and may be modified in non-backwards compatible
    |      ways or removed without deprecation or other notices.
  """)
class PostProcessGenotypeConcordance
( @arg(flag="c", doc="The VCF output of the genotype concordance tool.") val concordanceVcf: PathToVcf,
  @arg(flag="v", doc="The called VCF that went into genotype concordance.") val callsVcf: PathToVcf,
  @arg(flag="o", doc="Output tab-delimited text file.") val output: FilePath
) extends FgBioTool with LazyLogging {
  Io.assertReadable(concordanceVcf)
  Io.assertReadable(callsVcf)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val concIn = new VCFFileReader(concordanceVcf.toFile, false)
    val callIn = new VCFFileReader(callsVcf.toFile, false)
    val iterators = Seq(concIn, callIn).map(v => v.toIterator)
    val iterator  = JointVariantContextIterator(iterators, concIn.getFileHeader.getSequenceDictionary)

    val metrics = iterator.flatMap {
      case Seq(Some(conc), callOption) =>
        val callGt     = callOption.map(_.getGenotype(0)).getOrElse(conc.getGenotype("call"))
        val calledGtAd = Option(callGt.getAD)
        val concState  = conc.getAttributeAsStringList("CONC_ST", "").mkString(",")

        Some(VariantInfo(
          chrom        = conc.getContig,
          chr          = conc.getContig.replace("chr", ""),
          pos          = conc.getStart,
          ref          = conc.getReference.getBaseString,
          alt          = conc.getAlternateAlleles.mkString(","),
          mutation     = conc.getReference.getBaseString + ">" + conc.getAlternateAlleles.mkString(","),
          variant_type = if (conc.isSNP) "snp" else if (conc.isMNP) "mnp" else "indel",
          qual         = callOption.map(_.getPhredScaledQual),
          filter       = callOption.map(_.getFilters.mkString(",")),
          conc_state   = concState,
          conc_desc    = concStateToDescription(concState),
          expected_af  = Option(conc.getGenotype("truth").getExtendedAttribute("AF")).map(_.asInstanceOf[String].toDouble),
          genotype     = callGt.getGenotypeString(true),
          ref_obs      = calledGtAd.map(_(0)),
          alt_obs      = calledGtAd.map(_(1)),
          af           = calledGtAd.map(ad => ad(1) / ad.sum.toDouble)
        ))
      case Seq(None, Some(call)) =>
        None
    }.toSeq

    Metric.write[VariantInfo](metrics, output)
  }

  private def concStateToDescription(state: String): String = state match {
    case "FP,TN" => "FP"
    case "TN,FN" => "FN"
    case "TP"    => "TP"
    case "TP,FP" => "TP"
    case "TP,TN" => "TP"
    case _       => "its_complicated"
  }
}
