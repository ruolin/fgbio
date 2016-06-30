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

import com.fulcrumgenomics.util.NumericTypes._
import htsjdk.samtools.SAMUtils
import net.jafama.FastMath._

object NumericTypes {
  /** A numeric type representing a `-10 * log(10)` scaled probability of error. */
  type PhredScore = Byte

  object PhredScore {
    val MinValue: PhredScore = 2.toByte

    /** The maximum Phred Value. */
    val MaxValue: PhredScore = SAMUtils.MAX_PHRED_SCORE.toByte

    /** The maximum values as a LogDouble. */
    private val MaxValueAsLogDouble: LogDouble = LogDouble.fromPhredScore(MaxValue)

    /** The precision to use when converting to an Integer. */
    val Precision = 0.001

    /** The value of log(10.0). */
    private val Ln10: LogDouble = log(10.0)

    /** Caps the scores to the range MinValue..MaxValue. */
    def cap(score: PhredScore): PhredScore = {
      if (score < MinValue) MinValue
      else if (score > MaxValue) MaxValue
      else score
    }

    /** Converts a probability (log-space) to a phred score. */
    def fromLogProbability(lnProbError: LogDouble): PhredScore = {
      if (lnProbError < MaxValueAsLogDouble) MaxValue
      else Math.floor(-10.0 * (lnProbError/ Ln10) + Precision).toByte
    }

    /** Converts a linear (i.e. non-log) probability to a phred score. */
    def fromLinearProbability(p: Double): PhredScore = Math.floor(-10.0 * log10(p)).toByte

    /** Converts a numeric phred score to a one-byte character with a +33 offset. */
    def toFastq(score: PhredScore): Byte = SAMUtils.phredToFastq(score).toByte
  }

  /** A type representing a natural log value stored in a Double. */
  type LogDouble = Double

  object LogDouble {
    val LnZero: LogDouble = Double.NegativeInfinity
    val LnOne:   LogDouble = 0.0
    val LnTwo:   LogDouble = toLogDouble(2)
    val LnThree: LogDouble = toLogDouble(3)
    val LnFour:  LogDouble = toLogDouble(4)
    val LnFive:  LogDouble = toLogDouble(5)
    val LnSix:   LogDouble = toLogDouble(6)
    val LnSeven: LogDouble = toLogDouble(7)
    val LnEight: LogDouble = toLogDouble(8)
    val LnNine:  LogDouble = toLogDouble(9)
    val LnTen:   LogDouble = toLogDouble(10)

    /** Converts the double to log-space. */
    def toLogDouble(value: Double): LogDouble = {
      if (value < 0) throw new IllegalArgumentException("Cannot use LogDouble to store values less than zero.")
      log(value)
    }

    /** Takes a phred score and converts it into the natural log of the probability of error. */
    def fromPhredScore(s: PhredScore): LogDouble = LnTen * s / -10

    /** Takes a phred score and converts it into the natural log of the probability of error. */
    def fromPhredScore(s: Int): LogDouble = fromPhredScore(s.toByte)

    /** Calculates the log of the mean of the non-log values. */
    def mean(values: Array[LogDouble]): LogDouble = sum(values) - log(values.length)

    /** Implementation of addition. */
    def add(a: LogDouble, b: LogDouble): LogDouble = {
      if (a.isNegInfinity) b
      else if (b.isNegInfinity) a
      else if (b < a) add(b, a)
      else a + log1p(exp(b - a)) // using log1p, otherwise it would be `a + log(1.0 + exp(b - a))`
    }

    /** Implementation of subtraction. */
    def sub(a: LogDouble, b: LogDouble): LogDouble = {
      if (b.isNegInfinity) a
      else if (a == b) Double.NegativeInfinity
      else if (a < b) throw new IllegalArgumentException("Subtraction will be less than zero.")
      else a + log1p(0.0 - exp(b - a)) // using log1p, otherwise it would be `a + log(1.0 - exp(b - a))`
    }

    /** 1.0 - this; useful for probabilities. */
    def oneMinus(a: LogDouble): LogDouble = {
      if (LnOne < a) LnZero
      else if (a < LnZero) LnOne
      else sub(LogDouble.LnOne, a)
    }

    /** Calculated the sum of one or more LogDoubles. */
    def sum(values: Array[LogDouble]): LogDouble = {
      val (minValue, minValueIndex) = MathUtil.minWithIndex(values)
      if (minValue.isNegInfinity) minValue
      else {
        var sum: LogDouble = 0.0
        var i = 0
        while (i < values.length) {
          if (minValueIndex != i) sum += exp(values(i) - minValue)
          i += 1
        }
        minValue + log1p(sum)
      }
    }
  }
}

