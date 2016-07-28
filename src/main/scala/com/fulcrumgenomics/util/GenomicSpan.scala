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
 *
 */

package com.fulcrumgenomics.util

import htsjdk.samtools.util.{CoordMath, Locatable}

/** A simple trait that extends [[htsjdk.samtools.util.Locatable]] but provides a few easier to use
  * methods for getting the location.
  */
trait GenomicSpan extends Locatable {
  def contig: String
  def start: Int
  def end: Int
  def length: Int = end - start + 1

  /** True if the two genomic spans overlap, false otherwise. */
  def overlaps(other: GenomicSpan): Boolean = {
    if (this.contig != other.contig) false
    else CoordMath.overlaps(this.start, this.end, other.start, other.end)
  }

  /** True if this genomic span encloses the other genomic span, false otherwise. */
  def encloses(other: GenomicSpan): Boolean = {
    CoordMath.encloses(this.start, this.end, other.start, other.end)
  }

  def getContig: String = contig
  def getStart: Int = start
  def getEnd: Int = end
}
