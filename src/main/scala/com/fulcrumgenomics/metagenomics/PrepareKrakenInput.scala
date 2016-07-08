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

package com.fulcrumgenomics.metagenomics

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.fastq.{FastqRecord, FastqSource, FastqWriter}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import dagr.commons.util.LazyLogging
import dagr.sopt.cmdline.ValidationException
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMUtils
import htsjdk.samtools.util.TrimmingUtil


@clp(group=ClpGroups.Metagenomics, description=
"""
  |Takes one or more fastq files, and prepares a single fasta for kraken[1]. Individual reads
  |are quality trimmed using the phred/bwa method, and then any remaining low quality bases
  |are masked to Ns.  Read 1s and read 2s are the joined with a block of Ns between them
  |(of length kmerLength) so that the pair can be fed to kraken as a single read.
  |
  |Reads shorter than minReadLength after trimming are excluded from the output. For paired reads
  |if either read 1 or read 2 of a pair is shorter than minReadLength the pair is excluded.
  |
  |Input fastq files may be gzipped. Read 1 file(s) are required, read 2 file(s) are optional
  |but if provided must be equal in number to the read 1 files.
  |
  |[1] Kraken: http://ccb.jhu.edu/software/kraken/
""")
class PrepareKrakenInput
(
  @arg(flag="1", doc="One or more read1 fastq files.")         val read1: Seq[PathToFastq],
  @arg(flag="2", doc="One or more read2 fastq files.")         val read2: Seq[PathToFastq] = Seq.empty,
  @arg(flag="o", doc="Output (singular) fastq file to write.") val output: PathToFasta = Io.StdOut,
  @arg(flag="q", doc="Trim/mask bases below quality q.")       val minQuality: Int = 20,
  @arg(flag="k", doc="Kmer size of kraken database.")          val kmerSize: Int = 31,
  @arg(flag="l", doc="Discard reads shorter than this length.") val minReadLength: Int = 75
) extends FgBioTool with LazyLogging {
  if (read2.nonEmpty && read1.size != read2.size) {
    throw new ValidationException("If read2 fastqs are supplied, must be equal in number to read1 fastqs.")
  }

  Io.assertReadable(read1 ++ read2)
  Io.assertCanWriteFile(output)

  private val q2Phred    = SAMUtils.phredToFastq(2).toByte
  private val baseSpacer = "N" * kmerSize
  private val qualSpacer = q2Phred.toChar.toString * kmerSize
  private val minQualityPhred33 = SAMUtils.phredToFastq(minQuality).toByte

  override def execute(): Unit = {
    val fq1 = FastqSource(read1.map(Io.readLines).foldLeft(Iterator.empty.asInstanceOf[Iterator[String]])(_ ++ _))
    val fq2Maybe = if (read2.isEmpty) None else Some(FastqSource(read2.map(Io.readLines).foldLeft(Iterator.empty.asInstanceOf[Iterator[String]])(_ ++ _)))
    val out = Io.toWriter(output)
    val progress = new ProgressLogger(logger)

    // Utility method to write a fasta record to an output writer
    def writeFasta(name: String, bases: String): Unit = {
      out.write('>')
      out.write(name)
      out.newLine()
      out.write(bases)
      out.newLine()
    }

    while (fq1.hasNext) {
      val r1 = fq1.next()
      val r1Bases = trimAndMask(r1)

      fq2Maybe match {
        case None      =>
          if (r1Bases.length >= minReadLength) writeFasta(r1.name, r1Bases)
        case Some(fq2) =>
          val r2 = fq2.next()
          if (r1.name != r2.name) fail(s"Fastq files out of sync. Read one: ${r1.name} arrived with read two: ${r2.name}")

          val r2Bases = trimAndMask(r2)
          if (r1Bases.length >= minReadLength && r2Bases.length >= minReadLength) writeFasta(r1.name, r1Bases + baseSpacer + r2Bases)
      }

      progress.record()
    }

    out.close()
    fq1.safelyClose()
    fq2Maybe.foreach(_.safelyClose())
  }

  /** Trims the end of the read based on quality and then masks any remaining low quality bases to Ns. */
  def trimAndMask(fq: FastqRecord): String = {
    val trimPoint = TrimmingUtil.findQualityTrimPoint(SAMUtils.fastqToPhred(fq.quals), minQuality)
    val bases = fq.bases.getBytes()
    val quals = fq.quals.getBytes()

    var i = 0
    while (i<trimPoint) {
      if (quals(i) < minQualityPhred33) bases(i) = 'N'
      i += 1
    }

    new String(bases, 0, trimPoint)
  }
}
