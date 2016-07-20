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

import com.fulcrumgenomics.FgBioDef._
import htsjdk.samtools.util.{CloseableIterator, IntervalList}
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

object MultiIntervalVcfIterator {
  /** Makes an iterator from a VCFFileReader and IntervalList. */
  def apply(in: VCFFileReader, intervals: IntervalList) = new MultiIntervalVcfIterator(in, intervals)
}

class MultiIntervalVcfIterator private (val in: VCFFileReader, intervals: IntervalList) extends Iterator[VariantContext] {
  val intervalIterator = intervals.uniqued(false).iterator()
  private var vcfIterator: Option[CloseableIterator[VariantContext]] = None
  advance()

  /** If the current iterator is exhausted, move to the next interval. */
  private def advance() {
    // Discard any exhausted iterator
    if (this.vcfIterator.isDefined && !this.vcfIterator.get.hasNext) {
      this.vcfIterator.foreach(_.safelyClose())
      this.vcfIterator = None
    }

    // Find the next interval that has some variants in it and store that iterator
    while (vcfIterator.isEmpty && this.intervalIterator.hasNext) {
      val interval = this.intervalIterator.next
      val iterator = this.in.query(interval.getContig, interval.getStart, interval.getEnd)
      if (iterator.hasNext) this.vcfIterator = Some(iterator)
    }
  }

  /** Fetches the next VariantContext. */
  override def next(): VariantContext = this.vcfIterator match {
    case None    => throw new NoSuchElementException
    case Some(i) => yieldAndThen(i.next()) { advance() }
  }

  override def hasNext: Boolean = this.vcfIterator.isDefined
}
