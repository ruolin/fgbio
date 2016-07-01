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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._

object MathUtil {
  def mean(bytes: Array[Byte]): Byte = {
    var sum: Int = 0
    var i: Int = 0
    while (i < bytes.length) {
      sum += bytes(i)
      i += 1
    }

    (sum / bytes.length).toByte
  }

  /** Finds the minimum element in the array and returns it with its index. */
  def minWithIndex(vs: Array[Double], excludeNegInfinity: Boolean = true): (Double,Int) = {
    var min     = Double.MaxValue
    var minIndex = -1
    val len = vs.length
    var idx = 0
    while (idx < len) {
      val v = vs(idx)
      if (v < min && (!excludeNegInfinity || Double.NegativeInfinity != v)) {
        min = v
        minIndex = idx
      }
      idx += 1
    }

    (min, minIndex)
  }


  /** Finds the maximum element in the array and returns it with its index. */
  def maxWithIndex(vs: Array[Double]): (Double,Int) = {
    var max      = Double.MinValue
    var maxIndex = -1
    val len = vs.length
    var idx = 0
    while (idx < len) {
      val v = vs(idx)
      if (v > max) {
        max = v
        maxIndex = idx
      }
      idx += 1
    }

    (max, maxIndex)
  }
}
