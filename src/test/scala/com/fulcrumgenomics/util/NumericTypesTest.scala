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

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.NumericTypes._

import scala.math._

/**
  * Tests for LogProbability and PhredValue.
  */
class NumericTypesTest extends UnitSpec {
  
  "PhredScore" should "convert probabilities to phred scores" in {
    PhredScore.fromLogProbability(LogProbability.toLogProbability(0.0)) shouldBe PhredScore.MaxValue
    PhredScore.fromLogProbability(LogProbability.toLogProbability(0.1)) shouldBe 10
    PhredScore.fromLogProbability(LogProbability.toLogProbability(0.5)) shouldBe 3
    PhredScore.fromLogProbability(LogProbability.toLogProbability(1.0)) shouldBe 0
  }

  it should "convert phred scores to probabilities" in {
    PhredScore.toLogProbability(Byte.MaxValue) shouldBe -29.24283 +- PhredScore.Precision
    PhredScore.toLogProbability(10) shouldBe LogProbability.toLogProbability(0.1) +- PhredScore.Precision
    PhredScore.toLogProbability(3) shouldBe LogProbability.toLogProbability(0.5011872) +- PhredScore.Precision
    PhredScore.toLogProbability(0) shouldBe LogProbability.toLogProbability(1.0) +- PhredScore.Precision
  }

  it should "display phred scores as integers" in {
    PhredScore.fromLogProbability(PhredScore.toLogProbability(10)) shouldBe 10
    PhredScore.fromLogProbability(PhredScore.toLogProbability(PhredScore.MaxValue+1)) shouldBe PhredScore.MaxValue
    PhredScore.fromLogProbability(PhredScore.toLogProbability(9)) shouldBe 9
  }

  "LogProbability" should "convert doubles to log-space doubles" in {
    import LogProbability.toLogProbability
    toLogProbability(0)   shouldBe Double.NegativeInfinity
    toLogProbability(0.1) shouldBe -2.302585 +- 0.00001
    toLogProbability(0.5) shouldBe -0.6931472 +- 0.00001
    toLogProbability(1.0) shouldBe 0.0
  }

  it should "throw an exception if a negative value is to be converted to a log-space doubles" in {
    an[Exception] should be thrownBy LogProbability.toLogProbability(-1)
  }
  
  it should "or in log space" in {
    exp(LogProbability.or(10, 10)) shouldBe exp(10)*2 +- 0.00001
    exp(LogProbability.or(10, 20)) shouldBe exp(10)+exp(20) +- 0.00001
    exp(LogProbability.or(20, 10)) shouldBe exp(20)+exp(10) +- 0.00001
    exp(LogProbability.or(10, LnZero)) shouldBe exp(10)
    exp(LogProbability.or(LnZero, 10)) shouldBe exp(10)
    exp(LogProbability.or(10, LogProbability.toLogProbability(0))) shouldBe exp(10)
  }

  it should "notOther in log space" in {
    val Seq(q10, q20) = Seq(10, 20).map(LogProbability.fromPhredScore)
    LogProbability.notOther(10, 10) shouldBe LnZero
    LogProbability.notOther(q10, q10) shouldBe LnZero
    exp(LogProbability.notOther(q10, q20)) shouldBe (0.1-0.01) +- 0.00001
    an[IllegalArgumentException] should be thrownBy LogProbability.notOther(q20, q10)
    LogProbability.notOther(LnTen, LnZero) shouldBe LnTen
    an[IllegalArgumentException] should be thrownBy LogProbability.notOther(LogProbability.toLogProbability(1.0), LogProbability.toLogProbability(1.0000000000000004))
    an[IllegalArgumentException] should be thrownBy LogProbability.notOther(LogProbability.toLogProbability(1.0), LogProbability.toLogProbability(-0.0000000000000004))
  }

  it should "1 - probability in log space" in {
    import LogProbability.{not, toLogProbability}
    val Seq(q10, q20) = Seq(10, 20).map(LogProbability.fromPhredScore)
    exp(LogProbability.not(q10)) shouldBe 0.9 +- 0.00001
    exp(LogProbability.not(q20)) shouldBe 0.99 +- 0.00001
    exp(LogProbability.not(toLogProbability(0.90))) shouldBe 0.1  +- 0.00001
    exp(LogProbability.not(toLogProbability(0.99))) shouldBe 0.01 +- 0.00001
    exp(LogProbability.not(LnZero))  shouldBe 1.0  +- 0.00001

    exp(LogProbability.not(toLogProbability(1.0000000000000004))) shouldBe 0.0 +- 0.00001
    exp(LogProbability.not(toLogProbability(0.0000000000000004))) shouldBe 1.0 +- 0.00001
  }

  it should "and multiple values in log space" in {
    LogProbability.andAll(Array(log(10), log(10), log(10))) shouldBe log(1000) +- 0.00001
  }

  it should "or multiple values in log space" in {
    LogProbability.orAll(Array(log(10), log(10), log(10))) shouldBe log(30) +- 0.00001
  }
}
