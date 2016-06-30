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
  * Tests for LogDouble and PhredValue.
  */
class NumericTypesTest extends UnitSpec {
  "PhredScore" should "convert probabilities to phred scores" in {
    PhredScore.fromLogProbability(LogDouble.toLogDouble(0.0)) shouldBe PhredScore.MaxValue
    PhredScore.fromLogProbability(LogDouble.toLogDouble(0.1)) shouldBe 10
    PhredScore.fromLogProbability(LogDouble.toLogDouble(0.5)) shouldBe 3
    PhredScore.fromLogProbability(LogDouble.toLogDouble(1.0)) shouldBe 0
  }

//  it should "convert phred scores to probabilities" in {
//    Double.PositiveInfinity.fromPhredScore.linearValue shouldBe 0.0 +- Precision
//    10.0.fromPhredScore.linearValue shouldBe 0.1 +- Precision
//    3.0103.fromPhredScore.linearValue shouldBe 0.5 +- Precision
//    0.0.fromPhredScore.linearValue shouldBe 1.0
//  }

//  it should "display phred scores as integers" in {
//    10.0.fromPhredScore.toPhredScoreInt shouldBe 10
//    100.fromPhredScore.toPhredScoreInt shouldBe 100
//    (10-Precision).fromPhredScore.toPhredScoreInt shouldBe 10
//    (10-(Precision*10)).fromPhredScore.toPhredScoreInt shouldBe 9
//    9.0001.fromPhredScore.toPhredScoreInt shouldBe 9
//
//    10.0.fromPhredScore.toPhredScoreString shouldBe "10"
//    100.0.fromPhredScore.toPhredScoreString shouldBe "100"
//    (10-Precision).fromPhredScore.toPhredScoreString shouldBe "10"
//    (10-(Precision*10)).fromPhredScore.toPhredScoreString shouldBe "9"
//    (9.0001).fromPhredScore.toPhredScoreString shouldBe "9"
//  }

  "LogDouble" should "convert doubles to log-space doubles" in {
    import LogDouble.toLogDouble
    toLogDouble(0)   shouldBe Double.NegativeInfinity
    toLogDouble(0.1) shouldBe -2.302585 +- 0.00001
    toLogDouble(0.5) shouldBe -0.6931472 +- 0.00001
    toLogDouble(1.0) shouldBe 0.0
  }

  it should "throw an exception if a negative value is to be converted to a log-space doubles" in {
    an[Exception] should be thrownBy LogDouble.toLogDouble(-1)
  }

  it should "add doubles in log space" in {
    exp(LogDouble.add(10, 10)) shouldBe exp(10)*2 +- 0.00001
    exp(LogDouble.add(10, 20)) shouldBe exp(10)+exp(20) +- 0.00001
    exp(LogDouble.add(20, 10)) shouldBe exp(20)+exp(10) +- 0.00001
    exp(LogDouble.add(10, LogDouble.LnZero)) shouldBe exp(10)
    exp(LogDouble.add(LogDouble.LnZero, 10)) shouldBe exp(10)
    exp(LogDouble.add(10, LogDouble.toLogDouble(0))) shouldBe exp(10)
  }

  it should "subtract doubles in log space" in {
    val Seq(q10, q20) = Seq(10, 20).map(LogDouble.fromPhredScore)
    LogDouble.sub(10, 10) shouldBe LogDouble.LnZero
    LogDouble.sub(q10, q10) shouldBe LogDouble.LnZero
    exp(LogDouble.sub(q10, q20)) shouldBe (0.1-0.01) +- 0.00001
    an[IllegalArgumentException] should be thrownBy LogDouble.sub(q20, q10)
    LogDouble.sub(LogDouble.LnTen, LogDouble.LnZero) shouldBe LogDouble.LnTen
    an[IllegalArgumentException] should be thrownBy LogDouble.sub(LogDouble.toLogDouble(1.0), LogDouble.toLogDouble(1.0000000000000004))
    an[IllegalArgumentException] should be thrownBy LogDouble.sub(LogDouble.toLogDouble(1.0), LogDouble.toLogDouble(-0.0000000000000004))
  }

  it should "1 - doubles in log space" in {
    import LogDouble.{oneMinus, toLogDouble}
    val Seq(q10, q20) = Seq(10, 20).map(LogDouble.fromPhredScore)
    exp(oneMinus(q10)) shouldBe 0.9 +- 0.00001
    exp(oneMinus(q20)) shouldBe 0.99 +- 0.00001
    exp(oneMinus(toLogDouble(0.90))) shouldBe 0.1  +- 0.00001
    exp(oneMinus(toLogDouble(0.99))) shouldBe 0.01 +- 0.00001
    exp(oneMinus(LogDouble.LnZero))  shouldBe 1.0  +- 0.00001

    exp(oneMinus(toLogDouble(1.0000000000000004))) shouldBe 0.0 +- 0.00001
    exp(oneMinus(toLogDouble(0.0000000000000004))) shouldBe 1.0 +- 0.00001
  }

  it should "compute the mean of doubles in log space" in {
    exp(LogDouble.mean(Array(LogDouble.toLogDouble(0.1)))) shouldBe 0.1  +- 0.00001
    exp(LogDouble.mean(Array(LogDouble.toLogDouble(0.1), LogDouble.toLogDouble(0.1)))) shouldBe 0.1  +- 0.00001
    exp(LogDouble.mean(Array(LogDouble.toLogDouble(0.4), LogDouble.toLogDouble(0.2)))) shouldBe 0.3  +- 0.00001
  }

  it should "compare doubles in log space" in {
    LogDouble.toLogDouble(0.01) should be <  LogDouble.toLogDouble(0.1)
    LogDouble.toLogDouble(0.1)  should be >  LogDouble.toLogDouble(0.01)
    LogDouble.toLogDouble(0.1)  should be >= LogDouble.toLogDouble(0.1)
  }

  it should "addAll" in {
    LogDouble.sum(Array(log(10), log(10), log(10))) shouldBe log(30)
  }
}
