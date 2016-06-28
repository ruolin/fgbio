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
import com.fulcrumgenomics.util.{MolecularBarcode, ReadSegment, ReadStructure, Template}
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
    |read structures specified.  The tags to store the molecular barcodes will be associated to the molecular barcode
    |segment in the read structure based on the order specified.  If only one molecular barcode tag is given, then the
    |molecular barcodes will be concatenated and stored in that tag. Otherwise the number of molecular barcodes in the
    |read structure should match the number of tags given. Each end of a pair will store the molecular barcode(s) for
    |both ends.  Therefore, the tags should be specified only once.  Optionally, the read names can be annotated with
    |the molecular barcode(s) directly.  In this case, the read name will be formatted "<NAME>+<UMIs1><UMIs2>" where
    |"<UMIs1>" is the concatenation of read one's molecular barcodes.  Similarly for "<UMIs2>".
    |
    |Mapping information will be unchanged, as such, this tool should not be used on reads that have been mapped since
    |it will lead to an BAM with inconsistent records.
  """,
  group = ClpGroups.SamOrBam)
class ExtractUmisFromBam
( @arg(flag = "i", doc = "Input BAM file.")                                      val input: PathToBam,
  @arg(flag = "o", doc = "Output BAM file.")                                     val output: PathToBam,
  @arg(flag = "r", doc = "The read structure, one per read in a template.")      val readStructure: Seq[String],
  @arg(flag = "b", doc = "SAM tags in which to store the molecular barcodes (one-per segment).")
                                                                                 val molecularBarcodeTags: Seq[String] = Seq.empty,
  @arg(flag = "a", doc = "Annotate the read names with the molecular barcodes. See usage for more details.") val annotateReadNames: Boolean = true
) extends FgBioTool with LazyLogging {
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
      throw new ValidationException("Either a single SAM tag. or the same # of SAM tags as molecular barcodes in the read structure,  must be given.")
    }
    // ensure no duplicate tags
    if (molecularBarcodeTags.distinct.size != molecularBarcodeTags.length) {
      throw new ValidationException("Duplicate molcular barcode tags found: " +
        molecularBarcodeTags.groupBy(identity).collect { case (x, List(_, _, _*)) => x }.mkString(", "))
    }
  }

  // split them into to tags for segments in read 1 vs read 2
  private val (molecularBarcodeTags1, molecularBarcodeTags2) = {
    val numMolecularBarcodesRead1 = rs1.molecularBarcode.length
    (molecularBarcodeTags.subList(0, numMolecularBarcodesRead1), molecularBarcodeTags.subList(numMolecularBarcodesRead1, molecularBarcodeTags.length))
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
    iterator.sliding(2,2).foreach {
      case Seq(r1: SAMRecord, r2: SAMRecord) =>
        // verify everything is in order
        Seq(r1, r2).foreach(r => {
          if (!r.getReadPairedFlag) throw new IllegalStateException(s"Read ${r.getReadName} was not paired")
          if (!r1.getReadUnmappedFlag) throw new IllegalStateException(s"Read ${r.getReadName} was not mapped")
        })
        if (!r1.getFirstOfPairFlag) throw new IllegalStateException(s"Read ${r1.getReadName} was not the first end of a pair")
        if (!r2.getSecondOfPairFlag) throw new IllegalStateException(s"Read ${r2.getReadName} was not the second end of a pair")
        if (!r1.getReadName.equals(r2.getReadName)) throw new IllegalStateException(s"Read names did not match: '${r1.getReadName}' and '${r2.getReadName}'")
        // now do some work
        val bases1 = annotateRecord(record=r1, readStructure=rs1, molecularBarcodeTags=molecularBarcodeTags1)
        val bases2 = annotateRecord(record=r2, readStructure=rs2, molecularBarcodeTags=molecularBarcodeTags2)
        if (annotateReadNames) {
          // Developer Note: the delimiter must be an ascii character less than what is usually in the read names.  For
          // example "|" doesn't work.  I am not completely sure why.
          r1.setReadName(r1.getReadName + "+" + bases1 + bases2)
          r2.setReadName(r2.getReadName + "+" + bases1 + bases2)
        }
        assert(r1.getReadName.equals(r2.getReadName))
        molecularBarcodeTags1.foreach { tag => r2.setAttribute(tag, r1.getStringAttribute(tag)) }
        molecularBarcodeTags2.foreach { tag => r1.setAttribute(tag, r2.getStringAttribute(tag)) }
        out.addAlignment(r1)
        out.addAlignment(r2)
    }
    out.close()
    CloserUtil.close(iterator)
  }

  /**
    * Extracts bases associated with molecular barcodes and adds them as SAM tags.  The read's bases will only contain
    * template bases as defined in the given read structure.
    */
  private def annotateRecord(record: SAMRecord,
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
      case Seq(tag) => record.setAttribute(tag, molecularBarcodeBases.mkString)
      case _ =>
        assert(molecularBarcodeTags.length == molecularBarcodeBases.length)
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
