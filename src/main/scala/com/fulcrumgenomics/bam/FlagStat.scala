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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef.PathToBam
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.SafelyClosable
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Metric

/** Stores various metrics derived from the FLAG field in a SAM or BAM file.  The metric counts either reads that
  * pass-filter or reads that do not pass filter, but not both.
  *
  * See the [SAMTools flagstat documentation](http://www.htslib.org/doc/samtools.html)
  *
  * @param pf true if the metric is for pass-filter reads, false for reads that did not pass-filter.
  * @param reads the total number of reads/records.
  * @param mapped the number of mapped reads/records.
  * @param paired the number of reads/records that were paired.
  * @param mapped_in_pair the number of mapped paired reads/records with its mate mapped.
  * @param proper_paired the number of paired reads/records that were properly paired.
  * @param singleton the number of paired reads/records that were mapped but its mate was not mapped.
  * @param read_one the number of first of pair reads/records.
  * @param read_two the number of second of pair reads/records.
  * @param duplicate the number of duplicate reads/records.
  * @param inter_chrom the number of mapped paired reads/records whose mate maps to a different chromosome.
  * @param inter_chrom_mapq5 the number of mapped paired reads/records with mapping quality at least five whose mate maps to a different chromosome.
  * @param secondary the number of secondary reads/records.
  * @param supplementary the number of supplementary reads/records.
  */
case class FlagStatMetric
(pf: Boolean,
 var reads: Long = 0,
 var mapped: Long = 0,
 var paired: Long = 0,
 var mapped_in_pair: Long = 0,
 var proper_paired: Long = 0,
 var singleton: Long = 0,
 var read_one: Long = 0,
 var read_two: Long = 0,
 var duplicate: Long = 0,
 var inter_chrom: Long = 0,
 var inter_chrom_mapq5: Long = 0,
 var secondary: Long = 0,
 var supplementary: Long = 0
) extends Metric

@clp(group = ClpGroups.SamOrBam, description=
  """|Collects metrics from the FLAG field in a SAM or BAM file.
  """)
class FlagStat
( @arg(flag='i', doc="Input SAM or BAM file.") val input: PathToBam,
  @arg(flag='o', doc="The output file in which metrics will be stored.") val output: Option[PathToBam] = None
) extends FgBioTool with LazyLogging {

  override def execute(): Unit = {
    val metrics = Seq(true, false).map(pf => pf -> FlagStatMetric(pf=pf)).toMap
    val in      = SamSource(input)

    in.foreach { rec => FlagStat.updateMetric(metrics(rec.pf), rec) }
    in.safelyClose

    // Output to file, if given
    output.foreach { out => Metric.write(out, metrics.values.toSeq.sortBy(_.pf).reverse) } // NB: write pf reads first

    // Output the same format as SAMTools
    System.out.print(FlagStat.formatInSamtools(metrics))
  }
}

object FlagStat {
  /** Updates the metric using the given record. */
  private[bam] def updateMetric(metric: FlagStatMetric, rec: SamRecord): Unit = {
    require(metric.pf == rec.pf, "The metric PF value did not match the record's pf value")
    metric.reads += 1
    if (rec.secondary) metric.secondary += 1
    else if (rec.supplementary) metric.supplementary += 1
    else if (rec.paired) {
      metric.paired += 1
      if (rec.properlyPaired && rec.mapped) metric.proper_paired +=1
      if (rec.firstOfPair) metric.read_one += 1
      if (rec.secondOfPair) metric.read_two += 1
      if (rec.mapped && rec.mateUnmapped) metric.singleton += 1
      if (rec.mapped && rec.mapped) {
        metric.mapped_in_pair += 1
        if (rec.refIndex != rec.mateRefIndex) {
          metric.inter_chrom += 1
          if (rec.mapq >= 5) metric.inter_chrom_mapq5 += 1
        }
      }
    }
    if (rec.mapped) metric.mapped += 1
    if (rec.duplicate) metric.duplicate += 1
  }

  /** Creates a multi-line string equivalent to the output of SAMTools flagstat. */
  private[bam] def formatInSamtools(metrics: Map[Boolean, FlagStatMetric]): String = {
    require(metrics.size == 2, "Metrics did not contain both PF and non-PF metrics")

    /** Function to print the a metric values (ex. `read_one`) for both the PF and non-PF metrics. Prints a string with
      * the two values separated by a space-padded plus-sign (ex. `12 + 24`).*/
    def format_metric(func: FlagStatMetric => Long): String = f"${func(metrics(true))}%d + ${func(metrics(false))}%d"

    /** Function to create a string representing a percentage calculate from two metric values (ex. `mapped` divided by
      * `reads` for the percentage of mapped reads) for both the PF and non-PF metrics.  The values will be enclosed in
      * parentheses and a separated by " : " (ex. `(22.2 : 12.9)`. */
    def format_percent_string(numerator: FlagStatMetric => Long, denominator: FlagStatMetric => Long): String = {
      val percents = Seq(true, false).map { pf =>
        val m = metrics(pf)
        val d = denominator(m)
        if (d == 0) "N/A" else f"${100.0 * numerator(m) / d.toDouble}%.2f%%"
      }
      "(" + percents.mkString(" : ") + ")"
    }
    
    val builder = new StringBuilder

    /** Print the SAMTools flagstat output to standard output. */
    builder ++= format_metric(_.reads)             + " in total (QC-passed reads + QC-failed reads)" + "\n"
    builder ++= format_metric(_.secondary)         + " secondary" + "\n"
    builder ++= format_metric(_.supplementary)     + " supplementary" + "\n"
    builder ++= format_metric(_.duplicate)         + " duplicates" + "\n"
    builder ++= format_metric(_.reads)             + " mapped " + format_percent_string(_.mapped, _.reads) + "\n"
    builder ++= format_metric(_.paired)            + " paired in sequencing" + "\n"
    builder ++= format_metric(_.read_one)          + " read1" + "\n"
    builder ++= format_metric(_.read_two)          + " read2" + "\n"
    builder ++= format_metric(_.proper_paired)     + " properly paired " + format_percent_string(_.proper_paired, _.paired) + "\n"
    builder ++= format_metric(_.mapped_in_pair)    + " with itself and mate mapped" + "\n"
    builder ++= format_metric(_.singleton)         + " singletons " + format_percent_string(_.singleton, _.paired) + "\n"
    builder ++= format_metric(_.inter_chrom)       + " with mate mapped to a different chr" + "\n"
    builder ++= format_metric(_.inter_chrom_mapq5) + " with mate mapped to a different chr (mapQ>=5)" + "\n"
    
    builder.toString
  }
}


