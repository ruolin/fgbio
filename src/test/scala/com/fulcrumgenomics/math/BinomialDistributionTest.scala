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

import com.fulcrumgenomics.testing.UnitSpec
import org.apache.commons.math3.distribution.{BinomialDistribution => CommonsBinomial}

class BinomialDistributionTest extends UnitSpec {
  val binom = new BinomialDistribution

  "BinomialDistribution.coefficient" should "calculate the correct number of combinations given k=0 or k=n" in {
    binom.coefficient(100, 0)   shouldBe 1
    binom.coefficient(100, 100) shouldBe 1
  }

  it should "fail if n or k is less than 0" in {
    an[Exception] shouldBe thrownBy { binom.coefficient(100, -2) }
    an[Exception] shouldBe thrownBy { binom.coefficient(-2,  0) }
  }

  it should "fail if k > n" in {
    an[Exception] shouldBe thrownBy { binom.coefficient(n=100, k=101) }
  }

  Seq((20, 5, BigInt(15504)), (75, 25, BigInt("52588547141148893628"))).foreach { case (n, k, coeff) =>
    it should s"correctly calculate coefficient(n=$n, k=$k) = $coeff" in {
      binom.coefficient(n=n, k=k) shouldBe coeff
    }
  }

  it should "work for small values of n and k" in {
    binom.coefficient(20, 5) shouldBe BigInt("15504")
  }

  it should "work for moderate values of n and k" in {
    binom.coefficient(75, 25) shouldBe BigInt("52588547141148893628")
  }

  it should "work for very large values" in {
    binom.coefficient(200, 100) shouldBe BigInt("90548514656103281165404177077484163874504589675413336841320")
  }

  "BinomialCoefficient.cumulativeProbability" should "match commons math for small numbers" in {
    val commons = new CommonsBinomial(20, 0.05)
    Range.inclusive(0, 20).foreach { k =>
      binom.cumulativeProbability(k, n=20, p=0.05).toDouble shouldBe commons.cumulativeProbability(k) +- 1e-10
    }
  }

  it should "compute values where commons math runs out of precision" in {
    val commons = new CommonsBinomial(500, 0.01)
    commons.cumulativeProbability(250) shouldBe 1.0
    binom.cumulativeProbability(k=250, n=500, p=1.0) < BigDecimal(1) shouldBe true
  }
}
