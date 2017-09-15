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
 *
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, PathToBam, PathToFasta}
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.ProgressLogger
import enumeratum.EnumEntry
import com.fulcrumgenomics.FgBioDef.SafelyClosable
import com.fulcrumgenomics.bam.ReadEnd.{Both, First, Second}

import scala.collection.immutable.IndexedSeq


sealed trait ReadEnd extends EnumEntry
object ReadEnd extends FgBioEnum[ReadEnd] {
  def values: IndexedSeq[ReadEnd] = findValues
  case object First extends ReadEnd
  case object Second extends ReadEnd
  case object Both extends ReadEnd
}

@clp(group = ClpGroups.SamOrBam, description=
  """
    |Converts a paired end SAM/BAM into single-end data SAM/BAM.
    |
    |This tool will remove all pairing information.
  """)
class ConvertPairedEndToSingleEnd
( @arg(flag='i', doc="Input SAM or BAM file.") val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag='r', doc="Reference sequence fasta file.") val ref: Option[PathToFasta] = None,
  @arg(flag='e', doc="Which end to keep.") val end: ReadEnd = Both,
  @arg(flag='t', doc="The SAM tags that should be removed") val tags: Seq[String] = Seq("MC", "MQ", "Q2", "R2")
) extends FgBioTool with LazyLogging {
  override def execute(): Unit = {
    val in       = SamSource(input)
    val progress = ProgressLogger(logger)
    val out      = SamWriter(output, in.header, ref=ref)

    in.iterator
      .filter(keepEnd)
      .map(convert)
      .foreach { rec =>
        out += rec
        progress.record(rec)
      }

    out.close()
    in.safelyClose()
  }

  private def keepEnd(rec: SamRecord): Boolean = this.end match {
    case First  => rec.paired && rec.firstOfPair
    case Second => rec.paired && rec.secondOfPair
    case Both   => true
  }

  private def convert(rec: SamRecord): SamRecord = {
    // remove all pairing information
    if (rec.paired) rec.name = rec.name + (if (rec.firstOfPair) "/1" else "/2")
    rec.paired             = false
    rec.properlyPaired     = false
    rec.firstOfPair        = false
    rec.secondOfPair       = false
    rec.mateUnmapped       = false
    rec.mateRefIndex       = -1
    rec.matePositiveStrand = true
    rec.mateStart          = 0
    rec.insertSize         = 0

    // remove all mate tags
    this.tags.foreach { tag => rec(tag) = null }

    rec
  }
}