/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

import java.io.PrintStream
import java.security.SecureRandom

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Aligner
import com.fulcrumgenomics.alignment.Mode.Global
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.{DelimitedDataParser, LazyLogging, NumericCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Io

@clp(group=ClpGroups.Personal, description=
  """
    |Checks to see how random the tail of a cell barcode distribution is.
  """)
class CheckBarcodeRandomness
( @arg(flag='i', doc="Input file of barcodes and counts.") val input: FilePath,
  @arg(flag='o', doc="Output file prefix.") val output: PathPrefix,
  @arg(flag='t', doc="How many barcodes to take as good.") val take: Int,
  @arg(flag='d', doc="How many barcodes to drop before sampling.") val drop: Int,
  @arg(flag='s', doc="Sample size.") val sampleSize: Int = 1000,
  @arg(flag='T', doc="Number of threads to use.") val threads: Int = 10

) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  private val random  = new scala.util.Random(42)
  private val aligner = Aligner(matchScore=4, mismatchScore= -4, gapOpen= -3, gapExtend= -1, mode=Global)
  private val Bases = "ACGT".getBytes

  override def execute(): Unit = {
    val in = DelimitedDataParser(input, '\t')
    val cbcs = in.map { row => row[String]("cell_barcode") }.toIndexedSeq
    val random = new SecureRandom()

    val goodCbcs   = cbcs.take(this.take).map(_.getBytes)
    val badCbcs    = scala.util.Random.shuffle(cbcs.drop(take + drop)).take(sampleSize).map(_.getBytes)
    val randomCbcs = Range(0, sampleSize).map(i => randomer(cbcs.head.length))

    val badCounter    = new NumericCounter[Int]
    val randomCounter = new NumericCounter[Int]

    logger.info("Aligning 'bad' CBCs.")
    align(badCbcs, goodCbcs, badCounter)
    logger.info(s"Bad CBCs : mean=${badCounter.mean()}, stdev=${badCounter.stddev()}.")

    logger.info("Aligning randomers.")
    align(randomCbcs, goodCbcs, randomCounter)
    logger.info(f"Randomers: mean=${randomCounter.mean()}, stdev=${randomCounter.stddev()}.")

    val out = new PrintStream(output + ".txt")
    out.println("dataset\tscore\tcount")
    badCounter.foreach    { case (s, n) => out.println(s"bad\t$s\t$n") }
    randomCounter.foreach { case (s, n) => out.println(s"random\t$s\t$n") }
    out.close()
  }

  private def align(queries: Seq[Array[Byte]], targets: Seq[Array[Byte]], counter: NumericCounter[Int]): Unit = {
    queries.parWith(threads)
      .map { q => targets.map(t => aligner.align(q, t).score).max }
      .seq
      .foreach(score => counter.count(score))
  }

  /* Generates a random kmer of a given length. */
  private def randomer(len: Int): Array[Byte] = {
    val bytes = new Array[Byte](len)
    forloop(from=0, until=len) { i =>
      bytes(i) = Bases(this.random.nextInt(Bases.length))
    }

    bytes
  }
}


