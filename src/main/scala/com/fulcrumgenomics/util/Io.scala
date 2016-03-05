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
 */package com.fulcrumgenomics.util

import java.io.{OutputStream, InputStream}
import java.nio.file.{Files, Path}
import java.util.zip.{GZIPOutputStream, GZIPInputStream}

import dagr.commons.io.{PathUtil, IoUtil}

/**
  * Common IO Utility methods.
  */
object Io extends IoUtil {
  val BufferSize = 16 * 1024

  /** Adds the automatic handling of gzipped files when opening files for reading. */
  override def toInputStream(path: Path): InputStream = {
    PathUtil.extensionOf(path) match {
      case Some(".gz") => new GZIPInputStream(Files.newInputStream(path), BufferSize)
      case _           => super.toInputStream(path)
    }
  }

  /** Adds the automatic handling of gzipped files when opening files for writing. */
  override def toOutputStream(path: Path): OutputStream = {
    PathUtil.extensionOf(path) match {
      case Some(".gz") => new GZIPOutputStream(Files.newOutputStream(path), BufferSize)
      case _           => super.toOutputStream(path)
    }
  }
}
