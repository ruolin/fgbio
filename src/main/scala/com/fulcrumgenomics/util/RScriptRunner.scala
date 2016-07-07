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

import java.nio.file.{Files, Path, Paths}

import com.fulcrumgenomics.FgBioDef._

import scala.io.Source
import scala.sys.process.Process

object RScriptRunner {
  val R = "Rscript"

  /** Extracts an R script from a classpath resource to a file and runs it. */
  def run(resourcePath: String, args: Any*): Int = {
    val script = extract(resourcePath)
    yieldAndThen(run(script, args:_*)) { Files.delete(script) }
  }

  /** Runs an R script that is sitting on the fileystem. */
  def run(script: Path, args: Any*): Int = Process(command = R :: script.toString :: args.toList.map(_.toString))!

  /** Extracts a script as a resource from the classpath into a temporary file. */
  private[util] def extract(path: String): Path = {
    val in = getClass.getClassLoader.getResourceAsStream(path)
    if (in == null) throw new IllegalArgumentException(s"R Script not found in classpath: ${path}")

    val filePath = Files.createTempFile(Paths.get(path).getFileName.toString, ".R").toAbsolutePath
    val out = Io.toOutputStream(filePath)
    Io.copy(in, out)
    out.close()
    in.safelyClose()
    filePath
  }
}
