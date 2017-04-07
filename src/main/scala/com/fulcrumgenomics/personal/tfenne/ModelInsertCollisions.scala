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
import com.fulcrumgenomics.commons.util.{LazyLogging, NumericCounter, SimpleCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Metric, Rscript}
import org.apache.commons.math3.distribution.NormalDistribution

import scala.collection.mutable.ListBuffer

@clp(group=ClpGroups.Personal, description=
  """
    |Given a mean and standard deviation of an insert size distribution, model how frequently
    |we expect to see the same insert by chance.
  """)
class ModelInsertCollisions
( @arg(flag="m", doc="Mean of the insert size distribution.") val mean: Double,
  @arg(flag="s", doc="Standard deviation of the insert size distribution.") val sd: Double,
  @arg(flag="l", doc="The length of the chromosome to model.") chromLength: Int = 100000,
  @arg(flag="g", doc="Weight of one haploid copy of the genome in picograms.") genomeWeight: Double = 3,
  @arg(flag="d", doc="Amount of input DNA in nanograms.") inputDna: Seq[Double] = Seq(10.0),
  @arg(flag="o", doc="Prefix of output files to write.") val output: PathPrefix
) extends FgBioTool with LazyLogging {
  private val normal = new NormalDistribution(mean, sd)
  private val ScriptPath = "com/fulcrumgenomics/personal/tfenne/ModelInsertCollisions.R"

  /** Generates the next length that is greater than 0. */
  private def nextLength: Int = {
    val len = Math.round(normal.sample()).toInt
    if (len > 0) len else nextLength
  }

  override def execute(): Unit = {
    val allMetrics = ListBuffer[InsertCollisionMetric]()

    for (input <- inputDna) {
      val copies  = (input * 1000 / genomeWeight).toInt
      val counter = new SimpleCounter[String]
      logger.info(f"Simulating $copies%,d copies of the chromosome.")

      forloop (from=0, until=copies) { n =>
        var pos = nextLength

        while (pos < chromLength) {
          val start = pos
          val end   = start + nextLength - 1
          pos       = end + 1
          counter.count(start + ":" + end)
        }
      }

      // Generate a frequency of frequencies counter
      val freqOfFreqs = new NumericCounter[Long]
      counter.foreach { case (insert, count) => freqOfFreqs.count(count) }

      // Generate the metrics
      val total = freqOfFreqs.map(_._2).sum.toDouble
      var cumulative = 0L
      val metrics = freqOfFreqs.map { case (familySize, count) =>
        cumulative += count

        InsertCollisionMetric(
          input_dna            = input,
          genome_copies        = copies,
          mean_insert_size     = mean,
          sd_insert_size       = sd,
          family_size          = familySize.toInt,
          family_count         = count,
          fraction_of_families = count / total,
          cumulative_fraction  = cumulative / total
        )
      }

      allMetrics ++= metrics
    }

    val txtOutput = output.resolveSibling(output.getFileName + ".txt")
    val pdfOutput = output.resolveSibling(output.getFileName + ".pdf")

    Metric.write(txtOutput, allMetrics)
    Rscript.exec(ScriptPath, txtOutput.toString, pdfOutput.toString)
  }
}

case class InsertCollisionMetric
(
  input_dna: Double,
  genome_copies: Int,
  mean_insert_size: Double,
  sd_insert_size: Double,
  family_size: Int,
  family_count: Long,
  fraction_of_families: Double,
  cumulative_fraction: Double
) extends Metric
