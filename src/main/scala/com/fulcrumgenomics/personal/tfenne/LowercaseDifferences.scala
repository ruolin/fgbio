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

package com.fulcrumgenomics.personal.tfenne

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, SimpleCounter}
import dagr.sopt.{arg, clp}

@clp(group=ClpGroups.Personal, description="Lowercases bases that are different than the majority.")
class LowercaseDifferences(@arg(flag="i", doc="Input text file") val input: FilePath) extends FgBioTool {
  override def execute(): Unit = {
    val lines = Io.readLines(input).map(_.toUpperCase.toCharArray).toSeq
    val len = lines.map(_.length).max

    forloop(from = 0, until = len) { i =>
      val bases = lines.flatMap(l => if (i < l.length) Some(l(i)) else None)
      val counter = new SimpleCounter[Char]
      bases.foreach(counter.count)
      val (majorityBase, majorityCount) = counter.maxBy { case (base, count) => count }

      lines.foreach { line => if (line.length > i && line(i) != majorityBase) line(i) = line(i).toLower }
    }

    lines.foreach(line => println(new String(line)))
  }
}
