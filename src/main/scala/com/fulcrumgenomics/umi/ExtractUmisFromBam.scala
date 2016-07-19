/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{ProgressLogger, _}
import dagr.commons.CommonsDef.PathToBam
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt._
import dagr.sopt.cmdline.ValidationException
import htsjdk.samtools._
import htsjdk.samtools.util._

import scala.collection.JavaConversions._

@clp(description =
  """
    |Extracts unique molecular indexes from reads in a BAM file into tags.
    |
    |Currently only paired end read data are supported.
    |
    |Only template bases will be retained as read bases (stored in the SEQ field) as specified by the read structure.
    |
    |A read structure should be provided for each read of a template.  For example, paired end reads should have two
    |read structures specified.  The tags to store the molecular barcodes will be associated with the molecular barcode
    |segment(s) in the read structure based on the order specified.  If only one molecular barcode tag is given, then the
    |molecular barcodes will be concatenated and stored in that tag. Otherwise the number of molecular barcodes in the
    |read structure should match the number of tags given. Each end of a pair will contain the same molecular barcode
    |tags and values.
    |
    |Optionally, the read names can be annotated with the molecular barcode(s) directly.  In this case, the read name
    |will be formatted "<NAME>+<UMIs1><UMIs2>" where "<UMIs1>" is the concatenation of read one's molecular barcodes.
    |Similarly for "<UMIs2>".
    |
    |Mapping information will be unchanged, as such, this tool should not be used on reads that have been mapped since
    |it will lead to an BAM with inconsistent records.
    |
    |The read structure describes the structure of a given read as one or more read segments. A read segment describes
    |a contiguous stretch of bases of the same type (ex. template bases) of some length and some offset from the start
    |of the read.  The following segment types are supported:
    |  - T: template bases
    |  - B: sample barcode bases
    |  - M: molecular barcode bases
    |  - S: bases to ignore
    |An example would be "10B3M7S100T" which describes 120 bases, with the first ten bases being a sample barcode,
    |bases 11-13 being a molecular barcode, bases 14-20 ignored, and bases 21-120 being template bases.
  """,
  group = ClpGroups.SamOrBam)
class ExtractUmisFromBam
( @arg(flag = "i", doc = "Input BAM file.")                                      val input: PathToBam,
  @arg(flag = "o", doc = "Output BAM file.")                                     val output: PathToBam,
  @arg(flag = "r", doc = "The read structure, one per read in a template.")      val readStructure: Seq[String],
  @arg(flag = "b", doc = "SAM tags in which to store the molecular barcodes (one-per segment).")
                                                                                 val molecularBarcodeTags: Seq[String] = Seq.empty,
  @arg(flag = "a", doc = "Annotate the read names with the molecular barcodes. See usage for more details.") val annotateReadNames: Boolean = false
) extends FgBioTool with LazyLogging {

  val progress = new ProgressLogger(logger, verb="written", unit=5e6.toInt)
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  val (rs1, rs2) = readStructure match {
    case Seq(readStructure1, readStructure2) => (ReadStructure(readStructure1), ReadStructure(readStructure2))
    case _ => throw new ValidationException("Two read structures must be given (one per end of the template).")
  }

  // validate the read structure versus the molecular barcode tags
  {
    // create a read structure for the entire template
    val rs = ReadStructure(readStructure.mkString(""))
    // make sure each tag is of length 2
    molecularBarcodeTags.foreach(tag => if (tag.length != 2) throw new ValidationException("SAM tags must be of length two: " + tag))
    // ensure we either have one tag, or we have the same # of tags as molecular barcodes in the read structure.
    if (molecularBarcodeTags.size > 1 && rs.molecularBarcode.size != molecularBarcodeTags.size) {
      throw new ValidationException("Either a single SAM tag, or the same # of SAM tags as molecular barcodes in the read structure,  must be given.")
    }
  }

  // split them into to tags for segments in read 1 vs read 2
  private val (molecularBarcodeTags1, molecularBarcodeTags2) = molecularBarcodeTags match {
    case Seq(tag) => (Seq[String](tag), Seq[String](tag))
    case tags =>
      val numMolecularBarcodesRead1 = rs1.molecularBarcode.length
      (tags.subList(0, numMolecularBarcodesRead1).toSeq, tags.subList(numMolecularBarcodesRead1, tags.length).toSeq)
  }

  override def execute(): Unit = {
    assert(molecularBarcodeTags1.length == rs1.count {
      case m: MolecularBarcode => true
      case _ => false
    })
    assert(molecularBarcodeTags2.length == rs2.count {
      case m: MolecularBarcode => true
      case _ => false
    })
    val in: SamReader = SamReaderFactory.make.open(input.toFile)
    val iterator: SAMRecordIterator = in.iterator()
    val out: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader, true, output.toFile)
    iterator.grouped(2).foreach {
      case Seq(r1: SAMRecord, r2: SAMRecord) =>
        // verify everything is in order
        Seq(r1, r2).foreach(r => {
          if (!r.getReadPairedFlag) throw fail(s"Read ${r.getReadName} was not paired")
          if (!r.getReadUnmappedFlag) throw fail(s"Read ${r.getReadName} was not unmapped")
        })
        if (!r1.getFirstOfPairFlag) throw fail(s"Read ${r1.getReadName} was not the first end of a pair")
        if (!r2.getSecondOfPairFlag) throw fail(s"Read ${r2.getReadName} was not the second end of a pair")
        if (!r1.getReadName.equals(r2.getReadName)) throw fail(s"Read names did not match: '${r1.getReadName}' and '${r2.getReadName}'")
        // now do some work
        val bases1 = ExtractUmisFromBam.annotateRecord(record=r1, readStructure=rs1, molecularBarcodeTags=molecularBarcodeTags1)
        val bases2 = ExtractUmisFromBam.annotateRecord(record=r2, readStructure=rs2, molecularBarcodeTags=molecularBarcodeTags2)
        if (annotateReadNames) {
          // Developer Note: the delimiter must be an ascii character less than what is usually in the read names.  For
          // example "|" doesn't work.  I am not completely sure why.
          r1.setReadName(r1.getReadName + "+" + bases1 + bases2)
          r2.setReadName(r2.getReadName + "+" + bases1 + bases2)
        }
        assert(r1.getReadName.equals(r2.getReadName))

        // If we have duplicate tags, then concatenate them in the same order across the read pair.
        val tagAndValues = (molecularBarcodeTags1.map { tag => (tag, r1.getStringAttribute(tag)) }
          ++ molecularBarcodeTags2.map { tag => (tag, r2.getStringAttribute(tag)) }).toList
        tagAndValues.groupBy(_._1).foreach { case (tag, tuples) =>
            val attr = tuples.map(_._2).mkString(ExtractUmisFromBam.UmiDelimiter)
            r1.setAttribute(tag, attr)
            r2.setAttribute(tag, attr)
        }

        out.addAlignment(r1)
        out.addAlignment(r2)
        progress.record(r1, r2)
    }
    out.close()
    CloserUtil.close(iterator)
  }
}

object ExtractUmisFromBam {
  val UmiDelimiter = "-"

  /**
    * Extracts bases associated with molecular barcodes and adds them as SAM tags.  The read's bases will only contain
    * template bases as defined in the given read structure.
    */
  private[umi] def annotateRecord(record: SAMRecord,
                                  readStructure: ReadStructure,
                                  molecularBarcodeTags: Seq[String]): String = {
    val bases = record.getReadString
    val qualities = record.getBaseQualityString
    // get the bases associated with each segment
    val readStructureBases = readStructure.structureReadWithQualities(bases, qualities, strict = false)
    // get the molecular barcode segments
    val molecularBarcodeBases = readStructureBases.collect { case (b: String, q: String, m: MolecularBarcode) => b }
    // set the barcode tags
    molecularBarcodeTags match {
      case Seq(tag) => record.setAttribute(tag, molecularBarcodeBases.mkString(UmiDelimiter))
      case _ =>
        if (molecularBarcodeTags.length < molecularBarcodeBases.length) throw new IllegalStateException("Found fewer molecular barcode SAM tags than molecular barcodes in the read structure.")
        else if (molecularBarcodeTags.length > molecularBarcodeBases.length) throw new IllegalStateException("Found fewer molecular barcodes in the read structure than molecular barcode SAM tags.")
        molecularBarcodeTags.zip(molecularBarcodeBases).foreach { case (tag, b) => record.setAttribute(tag, b) }
    }
    // keep only template bases and qualities in the output read
    val basesAndQualities = readStructureBases.flatMap {
      case (b: String, q: String, s: ReadSegment) =>
        s match {
          case m: Template => Some((b, q))
          case _ => None
        }
    }
    record.setReadString(basesAndQualities.map(_._1).mkString)
    record.setBaseQualityString(basesAndQualities.map(_._2).mkString)
    // return the concatenation of the molecular barcodes in sequencing order
    molecularBarcodeBases.mkString
  }
}
