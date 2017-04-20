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

import java.util.Collections

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import dagr.commons.util.{BiMap, LazyLogging}
import dagr.sopt.{arg, clp}
import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory}
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.vcf.{VCFFileReader, VCFHeader, VCFHeaderLine}

@clp(group=ClpGroups.Personal, description=
  """
    |Converts a VCF between B37 and HG19 or vice versa. Not the most robust implementation, it
    |simply map the corresponding sequences in HG19 and B37 and then rewrites the VCF by:
    |  - Changing the sequence names in the records
    |  - Changing the contig lines in the header
    |  - Re-writing so that the records are in contig-order in case the ordering changed
    |
    |It does not do anyting clever around the different mitochondrial sequences.  And if your
    |reference (hg19 or B37) contains sequence names other than the expected set it will omit
    |those sequences!
  """)
class LiftVcfBetweenB37AndHg19
( @arg(flag="i", doc="Input VCF on either B37 or HG19.") val input: PathToVcf,
  @arg(flag="o", doc="Output VCF on on opposite build.") val output: PathToVcf,
  @arg(doc="Path to the HG19 fasta file.") val hg19: PathToFasta,
  @arg(doc="Path to the B37 fasta file.") val b37: PathToFasta
) extends FgBioTool with LazyLogging {

  Io.assertReadable(Seq(input, hg19, b37))
  Io.assertCanWriteFile(output)

  private val b37ToHg19 = Seq(
    "1"           ->  "chr1",
    "2"           ->  "chr2",
    "3"           ->  "chr3",
    "4"           ->  "chr4",
    "5"           ->  "chr5",
    "6"           ->  "chr6",
    "7"           ->  "chr7",
    "8"           ->  "chr8",
    "9"           ->  "chr9",
    "10"          ->  "chr10",
    "11"          ->  "chr11",
    "12"          ->  "chr12",
    "13"          ->  "chr13",
    "14"          ->  "chr14",
    "15"          ->  "chr15",
    "16"          ->  "chr16",
    "17"          ->  "chr17",
    "18"          ->  "chr18",
    "19"          ->  "chr19",
    "20"          ->  "chr20",
    "21"          ->  "chr21",
    "22"          ->  "chr22",
    "X"           ->  "chrX",
    "Y"           ->  "chrY",
    "MT"          ->  "chrM",
    "GL000207.1"  ->  "chr18_gl000207_random",
    "GL000226.1"  ->  "chrUn_gl000226",
    "GL000229.1"  ->  "chrUn_gl000229",
    "GL000231.1"  ->  "chrUn_gl000231",
    "GL000210.1"  ->  "chr21_gl000210_random",
    "GL000239.1"  ->  "chrUn_gl000239",
    "GL000235.1"  ->  "chrUn_gl000235",
    "GL000201.1"  ->  "chr9_gl000201_random",
    "GL000247.1"  ->  "chrUn_gl000247",
    "GL000245.1"  ->  "chrUn_gl000245",
    "GL000197.1"  ->  "chr8_gl000197_random",
    "GL000203.1"  ->  "chr17_gl000203_random",
    "GL000246.1"  ->  "chrUn_gl000246",
    "GL000249.1"  ->  "chrUn_gl000249",
    "GL000196.1"  ->  "chr8_gl000196_random",
    "GL000248.1"  ->  "chrUn_gl000248",
    "GL000244.1"  ->  "chrUn_gl000244",
    "GL000238.1"  ->  "chrUn_gl000238",
    "GL000202.1"  ->  "chr11_gl000202_random",
    "GL000234.1"  ->  "chrUn_gl000234",
    "GL000232.1"  ->  "chrUn_gl000232",
    "GL000206.1"  ->  "chr17_gl000206_random",
    "GL000240.1"  ->  "chrUn_gl000240",
    "GL000236.1"  ->  "chrUn_gl000236",
    "GL000241.1"  ->  "chrUn_gl000241",
    "GL000243.1"  ->  "chrUn_gl000243",
    "GL000242.1"  ->  "chrUn_gl000242",
    "GL000230.1"  ->  "chrUn_gl000230",
    "GL000237.1"  ->  "chrUn_gl000237",
    "GL000233.1"  ->  "chrUn_gl000233",
    "GL000204.1"  ->  "chr17_gl000204_random",
    "GL000198.1"  ->  "chr9_gl000198_random",
    "GL000208.1"  ->  "chr19_gl000208_random",
    "GL000191.1"  ->  "chr1_gl000191_random",
    "GL000227.1"  ->  "chrUn_gl000227",
    "GL000228.1"  ->  "chrUn_gl000228",
    "GL000214.1"  ->  "chrUn_gl000214",
    "GL000221.1"  ->  "chrUn_gl000221",
    "GL000209.1"  ->  "chr19_gl000209_random",
    "GL000218.1"  ->  "chrUn_gl000218",
    "GL000220.1"  ->  "chrUn_gl000220",
    "GL000213.1"  ->  "chrUn_gl000213",
    "GL000211.1"  ->  "chrUn_gl000211",
    "GL000199.1"  ->  "chr9_gl000199_random",
    "GL000217.1"  ->  "chrUn_gl000217",
    "GL000216.1"  ->  "chrUn_gl000216",
    "GL000215.1"  ->  "chrUn_gl000215",
    "GL000205.1"  ->  "chr17_gl000205_random",
    "GL000219.1"  ->  "chrUn_gl000219",
    "GL000224.1"  ->  "chrUn_gl000224",
    "GL000223.1"  ->  "chrUn_gl000223",
    "GL000195.1"  ->  "chr7_gl000195_random",
    "GL000212.1"  ->  "chrUn_gl000212",
    "GL000222.1"  ->  "chrUn_gl000222",
    "GL000200.1"  ->  "chr9_gl000200_random",
    "GL000193.1"  ->  "chr4_gl000193_random",
    "GL000194.1"  ->  "chr4_gl000194_random",
    "GL000225.1"  ->  "chrUn_gl000225",
    "GL000192.1"  ->  "chr1_gl000192_random"
  ).toMap

  private val hg19ToB37 = b37ToHg19.map { case (k,v) => (v,k) }

  override def execute(): Unit = {
    val hg19Ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(hg19)
    val b37Ref  = ReferenceSequenceFileFactory.getReferenceSequenceFile(b37)
    val in      = new VCFFileReader(input.toFile, true)

    val (mapping, sourceRef, targetRef) = in.iterator().next().getContig match {
      case chr if hg19ToB37.contains(chr) => (b37ToHg19, hg19Ref, b37Ref)
      case chr if b37ToHg19.contains(chr) => (hg19ToB37, b37Ref, hg19Ref)
      case chr => fail(s"Chromosome $chr not found in either HG19 or B37 sequence list.")
    }

    logger.info(s"Detected that input VCF is on ${sourceRef.getSequenceDictionary.getSequence(0).getAssembly}")

    val out = makeOutputWriter(in.getFileHeader, targetRef)
    val progress = new ProgressLogger(logger)

    // Loop through by target chromosome in order, pull out the variants and write to the target file
    targetRef.getSequenceDictionary.getSequences.foreach { targetChr =>
      val targetChrName = targetChr.getSequenceName
      val sourceChrName = mapping(targetChrName)
      val sourceChr     = sourceRef.getSequenceDictionary.getSequence(sourceChrName)
      val targetChrLength = targetChr.getSequenceLength

      val iterator      = in.query(sourceChr.getSequenceName, 1, sourceChr.getSequenceLength)
      iterator.filter(_.getEnd <= targetChrLength).foreach { sourceContext =>
        val builder = new VariantContextBuilder(sourceContext)
        builder.chr(targetChrName)
        val targetContext = builder.make()
        out.add(targetContext)
        progress.record(targetContext.getContig, targetContext.getStart)
      }

      iterator.close()
    }

    out.close()
  }

  /** Builds the output VCF writer with most of the input header, but switching sequence dictionary. */
  private def makeOutputWriter(in: VCFHeader, ref: ReferenceSequenceFile): VariantContextWriter = {
    val header = new VCFHeader(Collections.emptySet[VCFHeaderLine](), in.getSampleNamesInOrder)

    in.getFilterLines.foreach(header.addMetaDataLine)
    in.getFormatHeaderLines.foreach(header.addMetaDataLine)
    in.getInfoHeaderLines.foreach(header.addMetaDataLine)
    header.setSequenceDictionary(ref.getSequenceDictionary)

    val writer = new VariantContextWriterBuilder()
      .setOutputFile(output.toFile)
      .setOption(Options.INDEX_ON_THE_FLY)
      .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
      .setReferenceDictionary(ref.getSequenceDictionary)
      .build()

    writer.writeHeader(header)
    writer
  }
}
