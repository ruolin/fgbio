package com.fulcrumgenomics.rnaseq

import java.nio.file.Files

import com.fulcrumgenomics.commons.CommonsDef.{FilePath, PathPrefix}
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.commons.util.DelimitedDataParser
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}
import org.scalatest.OptionValues

import scala.io.Source

class AnnotateIntervalListTest extends UnitSpec with OptionValues {
  private def emtpyIntervalList(): IntervalList = {
    val header = new SAMFileHeader
    header.setSequenceDictionary(this.dict)
    new IntervalList(header)
  }

  "AnnotateIntervalList" should "annotate reads with annotations from and interval list" in {

  }

}
