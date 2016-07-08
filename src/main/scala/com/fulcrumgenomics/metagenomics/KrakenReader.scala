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

package com.fulcrumgenomics.metagenomics

import java.nio.file.Path

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.metagenomics.NcbiTaxonomy.TaxonId
import com.fulcrumgenomics.util.{Io, SimpleCounter}
import htsjdk.samtools.util.{StringUtil => HtsStringUtil}

/**
  * Represents a single result from a Kraken output file.
  *
  * @param queryname the queryname from the kraken file
  * @param taxonId   the assigned taxon ID from the kraken file
  * @param queryLength the query length from the taxon file
  * @param kmerCounts a map of taxon ID to kmers assigned to that taxon id, derived from the LCA
  *                   mapping information in the last column of the kraken output
  */
case class KrakenResult(queryname: String, taxonId: TaxonId, queryLength: Int, kmerCounts: Map[TaxonId,Int]) {
  def assigned: Boolean = this.taxonId != 0
}

/**
  * A reader for Kraken output files. If the path represents a regular file then the data
  * may be iterated more than once.  If the path represents a pipe or device then only a
  * single iterator can be performed.
  *
  * @param p the path to the kraken output file, optionally gzipped
  */
class KrakenReader(private val p: Path) extends Iterable[KrakenResult] {

  /** Generates an iterator over the kraken output. */
  override def iterator: Iterator[KrakenResult] = {
    Io.toSource(p).getLines().map { line =>
      val fields = new Array[String](5)
      HtsStringUtil.split(line, fields, '\t')

      val counter = new SimpleCounter[Int]
      fields(4).split(' ').filterNot(_.charAt(0) == 'A').map(_.split(':')).foreach(kv => counter.count(kv(0).toInt, kv(1).toInt))
      val counts = counter.iterator.map(ab => (ab._1, ab._2.toInt)).toMap

      KrakenResult(queryname=fields(1), taxonId=fields(2).toInt, queryLength=fields(3).toInt, kmerCounts=counts)
    }
  }
}
