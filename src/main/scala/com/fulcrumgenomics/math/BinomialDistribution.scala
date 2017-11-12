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

package com.fulcrumgenomics.math

import com.fulcrumgenomics.FgBioDef._

import scala.math.{BigDecimal, BigInt}

/**
  * A high-precision implementation of the math for binomial probabilities.  This implementation is several
  * times slower (though not orders of magnitude slower) than the implementation in Apache Commons Math, but
  * retains precision throughout the computation where the Commons implementation underflows.
  *
  * One implementation choice to be aware of when using this implementation is that each instance of the
  * class will calculate and cache factorials up to factorial(n) where n is the highest value of n
  * supplied to any call to [[probability()]] or [[cumulativeProbability()]].  Thus to ensure reasonable
  * performance it is recommended to create one instance of this class and reuse it for as may
  * calculations as is practical.
  */
class BinomialDistribution {
  // This array will get expanded with more factorials any time a higher n is queried
  private var factorials = Array(BigInt(0), BigInt(1))

  // Constant values for 0 and 1 in BigDecimal to avoid making them all the time
  private val Zero = BigDecimal(0)
  private val One  = BigDecimal(1)

  /** Get the factorial of n. Retrieves the value from a cached set of factorials. If the cache
    * does not contain factorial(n) yet, then all factorials up to factorial(n) are computed
    * and cached.
    *
    * @param n the number to compute the factorial of
    * @return the factorial value as a BigInt
    */
  private def factorial(n: Int): BigInt = {
    if (n > factorials.length - 1) {
      // Expand the array
      val oldLength = factorials.length
      val tmp = new Array[BigInt](n+1)
      System.arraycopy(factorials, 0, tmp, 0, oldLength)
      factorials = tmp

      // Fill in the new slots
      forloop (from=oldLength, until=factorials.length) { i =>
        factorials(i) = factorials(i-1) * BigInt(i)
      }
    }

    factorials(n)
  }

  /** Computes the binomial coefficient, i.e. the number of combinations of k items given a total of
    * n items, often phrased "n pick k".
    *
    * @param n The number of items being picked from
    * @param k The number of items being picked (must be 0 <= k <= n)
    * @return The number of combinations of k items from n
   */
  def coefficient(n: Int, k: Int): BigInt = {
    require(k <= n, s"k must be <= n. n=$n, k=$k")
    require(k >= 0, s"Can't compute coefficient for picking fewer than 0 values. n=$n, k=$k")

    if (k == 0 || k == n) 1
    else factorial(n) / (factorial(k) * factorial(n - k))
  }

  /**
    * Calculates the probability of exactly `k` successes in `n` trials given `p` probability.
    * Akin to `dbinom` in R.
    *
    * @param k The number of successful trials
    * @param n The number of trials
    * @param p The probability of the success in any one trial
    * @return The probability of exactly k successes in n trials of p probability
    */
  def probability(k: Int, n: Int, p: Double): BigDecimal = probability(k, n, BigDecimal(p))

  /**
    * Calculates the probability of exactly `k` successes in `n` trials given `p` probability.
    * Akin to `dbinom` in R.
    *
    * @param k The number of successful trials
    * @param n The number of trials
    * @param p The probability of the success in any one trial
    * @return The probability of exactly k successes in n trials of p probability
    */
  def probability(k: Int, n: Int, p: BigDecimal): BigDecimal = {
    require(p >= Zero && p <= One, s"Probability p must be between 0 and 1. p=$p.")
    val coeff = BigDecimal(coefficient(n, k))
    val p1 = p.pow(k)
    val p2 = (One-p).pow(n - k)
    coeff * p1 * p2
  }

  /**
    * Calculates the cumulative probability of `[0, k]` successes from n trials with `p` probability
    * of success in any individual trial.
    * Akin to `pbinom` in R.
    *
    * @param k the number of successes to compute the cumulative probability up to, inclusive
    * @param n the number of trials
    * @param p the probability of success in any single trial
    */
  def cumulativeProbability(k: Int, n: Int, p: Double): BigDecimal = cumulativeProbability(k, n, BigDecimal(p))

  /**
    * Calculates the cumulative probability of `[0, k]` successes from n trials with `p` probability
    * of success in any individual trial.
    * Akin to `pbinom` in R.
    *
    * @param k the number of successes to compute the cumulative probability up to, inclusive
    * @param n the number of trials
    * @param p the probability of success in any single trial
    */
  def cumulativeProbability(k: Int, n: Int, p: BigDecimal): BigDecimal = {
    var result = Zero
    forloop (from=0, until=k+1) { ki =>
      result += probability(ki, n, p)
    }

    result
  }
}
