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

import htsjdk.samtools.SAMUtils
import net.jafama.FastMath._

object NumericTypes {
  /** The value of toLogProbability(x) for [0, 1, ..., 10] */
  private val LnValues: Array[Double] = Array(Double.NegativeInfinity, 0.0) ++ Array.range(2, 11, 1).map(x => log(x))
  val LnZero   = LnValues(0)
  val LnOne    = LnValues(1)
  val LnTwo    = LnValues(2)
  val LnThree  = LnValues(3)
  val LnFour   = LnValues(4)
  val LnFiv    = LnValues(5)
  val LnSiz    = LnValues(6)
  val LnSeven  = LnValues(7)
  val LnEight  = LnValues(8)
  val LnNone   = LnValues(9)
  val LnTen    = LnValues(10)


  def lnValueDouble(value: Double): Double = log(value)

  def lnValueInt(value: Int): Double = {
    if (0 <= value && value < LnValues.length) LnValues(value)
    else log(value)
  }

  /** A numeric type representing a `-10 * log(10)` scaled probability of error. */
  type PhredScore = Byte

  object PhredScore {
    val MinValue: PhredScore = 2.toByte

    /** The maximum Phred Value. */
    val MaxValue: PhredScore = SAMUtils.MAX_PHRED_SCORE.toByte

    /** The maximum values as a LogDouble. */
    private val MaxValueAsLogDouble: LogProbability = LogProbability.fromPhredScore(MaxValue)

    /** The precision to use when converting to an Integer. */
    val Precision = 0.001

    /** The value of log(10.0). */
    private val Ln10: LogProbability = log(10.0)

    /** Caps the scores to the range MinValue..MaxValue. */
    def cap(score: PhredScore): PhredScore = {
      if (score < MinValue) MinValue
      else if (score > MaxValue) MaxValue
      else score
    }

    /** Converts a probability (log-space) to a phred score. */
    def fromLogProbability(lnProbError: LogProbability): PhredScore = {
      if (lnProbError < MaxValueAsLogDouble) MaxValue
      else Math.floor(-10.0 * (lnProbError/ Ln10) + Precision).toByte
    }

    /** Converts a linear (i.e. non-log) probability to a phred score. */
    def fromLinearProbability(p: Double): PhredScore = Math.floor(-10.0 * log10(p)).toByte

    /** Converts a numeric phred score to a one-byte character with a +33 offset. */
    def toFastq(score: PhredScore): Byte = SAMUtils.phredToFastq(score).toByte

    /** Takes a phred score and converts it into the natural log of the probability of error. */
    def toLogProbability(s: PhredScore): LogProbability = LnValues(10) * s / -10

    /** Takes a phred score and converts it into the natural log of the probability of error. */
    def toLogProbability(s: Int): LogProbability = {
      if (s > Byte.MaxValue) toLogProbability(Byte.MaxValue)
      else toLogProbability(s.toByte)
    }
  }

  /** A type representing a natural log value of a probability. */
  type LogProbability = Double

  object LogProbability {
    /** Converts the double to log-space. */
    def toLogProbability(value: Double): LogProbability = {
      if (value < 0) throw new IllegalArgumentException("Cannot use LogDouble to store values less than zero: " + value)
      log(value)
    }

    /** Takes a phred score and converts it into the natural log of the probability of error. */
    def fromPhredScore(s: PhredScore): LogProbability = LnValues(10) * s / -10

    /** Takes a phred score and converts it into the natural log of the probability of error. */
    def fromPhredScore(s: Int): LogProbability = fromPhredScore(s.toByte)

    /** Computes the probability of a or b, where a and b are independent events: Pr(A or B) = Pr(A) + Pr(B). */
    def or(a: LogProbability, b: LogProbability): LogProbability = {
      if (a.isNegInfinity) b
      else if (b.isNegInfinity) a
      else if (b < a) or(b, a)
      else a + log1p(exp(b - a)) // using log1p, otherwise it would be `a + log(1.0 + exp(b - a))`
    }

    /** Computes the probability of any of the given independent events occurring: Pr(AB..N) = Pr(A)+Pr(B)+...+Pr(N). */
    def orAll(values: Array[LogProbability]): LogProbability = {
      if (values.forall(_.isNegInfinity)) Double.NegativeInfinity
      else {
        val (minValue, minValueIndex) = MathUtil.minWithIndex(values)
        var sum: LogProbability = 0.0
        var i = 0
        while (i < values.length) {
          if (minValueIndex != i) sum += exp(values(i) - minValue)
          i += 1
        }
        minValue + log1p(sum)
      }
    }

    /** Computes the probability of a and b, where a and b are independent events: Pr(AB) = Pr(A)*Pr(B). */
    def and(a: LogProbability, b: LogProbability): LogProbability = a + b

    /** Computes the probability of the given independent events co-occurring: Pr(AB..N) = Pr(A)*Pr(B)*...*Pr(N). */
    def andAll(values: Array[Double]): LogProbability = values.sum

    /** Computes the probability Pr(A and not B) = Pr(A) - Pr(B). */
    def notOther(a: LogProbability, b: LogProbability): LogProbability = {
      if (b.isNegInfinity) a
      else if (a == b) Double.NegativeInfinity
      else if (a < b) throw new IllegalArgumentException("Subtraction will be less than zero.")
      else a + log1p(0.0 - exp(b - a)) // using log1p, otherwise it would be `a + log(1.0 - exp(b - a))`
    }

    /** Returns the probability of Pr(not A) = 1 - Pr(A).*/
    def not(a: LogProbability): LogProbability = {
      if (LnValues(1)  < a) LnValues(0)
      else if (a < LnValues(0) ) LnValues(1)
      else notOther(LnValues(1), a)
    }

    /** Normalizes a probability. */
    def normalizeByScalar(p: LogProbability, divideBy: Int): LogProbability = {
      if (0 <= divideBy && divideBy < LnValues.length) div(p, LnValues(divideBy))
      else div(p, toLogProbability(divideBy))
    }

    /** Normalizes a probability. */
    def normalizeByLogProbability(p: LogProbability, divideBy: LogProbability): LogProbability = div(p, divideBy)

    /** Divides a by b. */
    private def div(a: LogProbability, b: LogProbability): LogProbability = {
      if (b == Double.PositiveInfinity) throw new IllegalStateException("Trying to divide by 0: " + b)
      a - b
    }
  }
}

