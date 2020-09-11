/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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

package com.fulcrumgenomics.vcf.parsing.variant

import com.fulcrumgenomics.fasta.SequenceMetadata
import com.fulcrumgenomics.vcf.api.Allele.{NoCallAllele, SpannedAllele, SymbolicAllele}
import com.fulcrumgenomics.vcf.api.{Allele, AlleleSet, ArrayAttr, Genotype, Variant, VcfCount, VcfCountAndKindHeader, VcfFieldType, VcfFilterHeader, VcfFormatHeader, VcfHeader, VcfInfoHeader}
import fastparse.NoWhitespace._
import fastparse.{&, CharIn, CharPred, Index, P, _}
import com.fulcrumgenomics.vcf.parsing.util.ParseResult._
import com.fulcrumgenomics.vcf.parsing.util.ValueAndIndex._
import com.fulcrumgenomics.vcf.parsing.util.ValueAndIndex.CaptureResultAndIndex

import scala.annotation.tailrec
import scala.collection.immutable.ListMap
import scala.util.{Success, Try}


// TODO:
// - id must be unique across variants
// - coordinate sort order...
//

sealed trait BlahUtil {

  def field[_ :P]: P[SI] = (Index ~ CharPred(_ != '\t').rep(1).! ~ Index).si
  def tabThenField[_ :P]: P[SI] = "\t" ~/ field

}

/** Trait to parse the contig name and return the [[SequenceMetadata]] */
sealed trait SequenceMetadataParser extends BlahUtil {
  def header: VcfHeader

  /** Parses the contig name and returns the associated [[SequenceMetadata]], or fails if not found. */
  def parseSequenceMetadata[_: P]: P[SequenceMetadata] = (Index ~ field ~ Index).si.flatMap { name =>
    header.dict.get(name.value) match {
      case None       => name.fail(s"contig not found in header: `${name.value}`")
      case Some(data) => success(data)
    }
  }
}

sealed trait PositionParser {
  def parsePosition[_: P]: P[Int] = {
    (CharIn("1-9").! ~ CharsWhileIn("0-9").!).map { case (left, right) => (left + right).toInt }
  }
}

sealed trait IdentifierParser {

  private def parseId[_: P]: P[SI] = {
    (Index ~ CharsWhile(char => char != ';' && char != ' ').! ~ Index).si
  }

  def parseIdentifier[_: P]: P[Seq[SI]] = {
    (Index ~ ".".! ~ Index).si.map(si => Seq(si)) | parseId.rep(1)
  }
}

sealed trait RefAlleleParser extends BlahUtil {

  private val DnaBases: String = "ACGTNacgtn"

  def parseRefAllele[_: P](metadata: SequenceMetadata, pos: Int): P[String] = field.flatMap { refSI =>
    val ref: String = refSI.value
    val refEnd: Int = pos + ref.length
    if (refEnd > metadata.length) {
       refSI.fail(s"variant end was past the end of the reference: `$refEnd > ${metadata.length}`")
    }
    else {
      ref.indexWhere(base => !DnaBases.contains(base)) match {
        case -1    => success(ref)
        case index =>
          val prefix      = ref.take(index)
          val invalidBase = ref(index)
          val suffix      = ref.drop(index + 1)
          refSI.fail(s"invalidate ref base: `$prefix[$invalidBase]$suffix")
      }
    }
  }
}

sealed trait AltAlleleParser extends BlahUtil {

  def dna[_: P]: P[Allele] = CharIn("ACGTNacgtn").rep(1).!.map(Allele.apply)
  def spanendAllele[_: P]: P[Allele] = "*".!.map(_ => SpannedAllele)
  def missingAllele[_: P]: P[Allele] = ".".!.map(_ => NoCallAllele)
  def symbolicAllele[_: P]: P[Allele] = ("<" ~ CharPred(_ != '>').! ~ ">").map(SymbolicAllele.apply)
  // TODO: breakend

  def alts[_: P]: P[Seq[Allele]] = {
      dna.rep(1, sep=",") |
        spanendAllele.rep(exactly=1) |
        missingAllele.rep(exactly=1) |
        symbolicAllele.rep(exactly=1)
  }

  def parseAltAlleles[_: P](ref: String): P[AlleleSet] = alts.map { alleles =>
    AlleleSet(
      ref  = Allele(ref),
      alts = alleles
    )
  }
}

sealed trait QualParser {
  private def double[_: P]: P[Double] = {
    (CharIn("1-9") ~ CharsWhileIn("0-9").rep ~ ".".? ~ CharsWhileIn("0-9")).!.map(_.toDouble)
  }
  def parseQual[_: P]: P[Option[Double]] = ".".!.map(_ => None) | double.map(Some(_))
}

sealed trait FilterParser {
  def header: VcfHeader

  private val illegalFilterChars: String = "; \t\n"
  private def filterValue[_: P]: P[VcfFilterHeader] = {
    (Index ~ CharPred(!illegalFilterChars.contains(_)).rep(1).! ~ Index).si.flatMap { key =>
      header.filter.get(key.value) match {
        case Some(filter) => success(filter)
        case None         => key.fail(s"FILTER was not found in the header: `${key.value}`")
      }
    }
  }
  private def filterValues[_: P]: P[Seq[VcfFilterHeader]] = filterValue.rep(min=0, sep=":")
  private def missingValue[_: P]: P[Seq[VcfFilterHeader]] = ".".!.map(_ => Seq.empty)
  def parseFilter[_: P]: P[Seq[VcfFilterHeader]] = missingValue | filterValues
}

sealed trait Blah {
  protected case class Field[T](key: SI, field: VcfCountAndKindHeader, values: Seq[T] = Seq.empty)

}

sealed trait InfoParser extends Blah {
  def header: VcfHeader

  /** Parses the INFO key, making sure it is found in the VCF header */
  def infoKey[_: P]: P[Field[SI]] = {
    (Index ~ ("1000G".! | (CharIn("A-Za-z_") ~ CharIn("0-9A-Za-z_.").rep).!) ~ Index).si.flatMap { key =>
      this.header.info.get(key.value) match {
        case None       => key.fail(f"INFO key not found in the header: `${key.value}`")
        case Some(info) => success(Field[SI](key, info))
      }
    }
  }

  /** Parses the value for an INFO field. Must be `^([A-Za-z ][0-9A-Za-z .]*|1000G)$,` */
  def infoValue[_: P]: P[SI] = (Index ~ CharPred(c => c != ';' && c != '=').rep(1).! ~ Index).si

  /** Parses an INFO key and value (may be a list) */
  def infoKeyAndValues[_: P]: P[Field[SI]] = {
    (infoKey ~ ("=" ~ infoValue.rep(min=1, sep=",")).?).map {
      case (data, None)         => data
      case (data, Some(values)) => data.copy(values=values)
    }
  }

  def infos[_: P]: P[Seq[Field[SI]]] = {
     (".".!.map(_ => Seq.empty) | infoKeyAndValues.rep(sep=";"))
  }

  private def parseMissing[_: P](missing: Field[SI], infoValues: Seq[Field[SI]]): P[ListMap[String, Any]] = {
    if (infoValues.length != 1) {
      missing.key.fail("found `.` with other INFO fields")
    }
    else if (missing.values.nonEmpty) {
      missing.values.head.fail("INFO field `.` cannot have a value")
    }
    else {
      success(ListMap.empty)
    }
  }

  private def parseDuplicateKeys[_: P](infoValues: Seq[Field[SI]]): P[Seq[Field[SI]]] = {
    infoValues.groupBy(_.key.value).find(_._2.length > 1) match {
      case None                => success(infoValues)
      case Some((key, values)) => values.head.key.fail(f"found INFO key `$key` ${values.length} times")
    }
  }

  private def parseCount[_: P](infoValues: Seq[Field[SI]], alleles: AlleleSet): P[Seq[Field[SI]]] = {
    infoValues.find { case Field(_, info, values) =>
      val countOption: Option[Int] = info.count match {
        case VcfCount.Fixed(count)    => Some(count)
        case VcfCount.OnePerAllele    => Some(alleles.size)
        case VcfCount.OnePerAltAllele => Some(alleles.size - 1)
        case VcfCount.Unknown         => None
        case VcfCount.OnePerGenotype  => throw new Exception(s"bug: found one-per-genotype in INFO.${info.id}")
      }
      countOption.exists(count => count != values.length)
    } match {
      case None => success(infoValues)
      case Some(Field(key, info, values)) =>
        values.headOption.getOrElse(key).fail(
          s"INFO key `${key.value}`: expected ${info.count} values but found ${values.length} values"
        )
    }
  }

  sealed trait ParseResult {
    def isSuccess: Boolean
    def isFailure: Boolean = !isSuccess
    def success: ParseSuccess
    def failure: ParseFailure
  }
  private case class ParseSuccess(info: VcfCountAndKindHeader, results: Seq[Any]) extends ParseResult {
    override def isSuccess: Boolean = true
    def success: ParseSuccess = this
    def failure: ParseFailure = throw new NoSuchElementException
  }
  private case class ParseFailure(info: VcfCountAndKindHeader, value: SI, ex: Exception) extends ParseResult {
    override def isSuccess: Boolean = false
    def success: ParseSuccess = throw new NoSuchElementException
    def failure: ParseFailure = this

  }

  @tailrec
  private def parseRecursive(values: Seq[SI], curSuccess: ParseSuccess): ParseResult = {
      values match {
        case Seq()                  => curSuccess
        case Seq(value, tail @ _*) =>
          val parseResult = try {
            curSuccess.copy(
              results = curSuccess.results :+ curSuccess.info.kind.parse(value.value)
            )
          } catch {
            case ex: Exception => ParseFailure(info=curSuccess.info, value=value, ex=ex)
          }
          // Developer note: this is to make the function tail recursive
          parseResult match {
            case failure: ParseFailure => failure
            case success: ParseSuccess => parseRecursive(values=tail, curSuccess=success)
          }
      }
  }

  private def parseValues(info: VcfCountAndKindHeader, values: Seq[SI]): ParseResult = {
    val startSuccess = ParseSuccess(info=info, results=Seq.empty)
    values match {
      case Seq(value) if value.value == "." =>
        parseRecursive(values=Seq.empty, curSuccess=startSuccess)
      case _         =>
        parseRecursive(values=values, curSuccess=startSuccess)
    }
  }

  @tailrec
  private def parseKindRecursive(infoValues: Seq[Field[SI]], pastResults: Seq[ParseResult]): Seq[ParseResult] = {
    infoValues match {
      case Seq()                => pastResults
      case Seq(data, tail @ _*) =>
        parseValues(info=data.field, values=data.values) match {
          case failure: ParseFailure => Seq(failure)
          case success: ParseSuccess => parseKindRecursive(infoValues=tail, pastResults=pastResults :+ success)
        }
    }
  }

  private def parseKind[_: P](infoValues: Seq[Field[SI]]): P[Seq[(String,Any)]] = {
    parseKindRecursive(infoValues=infoValues, pastResults=Seq.empty) match {
      case Seq(ParseFailure(info, value, ex)) =>
        value.fail(s"could not parse `${value.value}` as type `${info.kind}`: ${ex.getMessage}")
      case results =>
        success(results.map(_.success).map {
          case ParseSuccess(info, Seq(result)) => info.id -> result
          case ParseSuccess(info, results)     => info.id -> results
        })
    }
  }

  def parseInfos[_: P](alleles: AlleleSet): P[ListMap[String, Any]] = infos.flatMap { infoValues=>
    infoValues.find(_.key.value == ".") match {
      case Some(missing) => parseMissing(missing=missing, infoValues=infoValues)  // empty INFO
      case None          => // non-empty INFOs
        parseDuplicateKeys(infoValues=infoValues)
          .flatMap { _ => parseCount(infoValues=infoValues, alleles=alleles) }
          .flatMap { _ => parseKind(infoValues=infoValues) }
          .map { values =>
            val builder = ListMap.newBuilder[String, Any]
            builder ++= values
            builder.result()
          }
    }
  }
}

sealed trait FormatKeysParser  {
  def header: VcfHeader

  def formatKey[_: P]: P[VcfFormatHeader] = (Index ~ (CharIn("A-Za-z_") ~ CharIn("0-9A-Za-z_.").rep).! ~ Index).si.flatMap {
    key: SI =>
      this.header.format.get(key.value) match {
        case None         => key.fail(f"FORMAT key not found in the header: `${key.value}`")
        case Some(format) => success(format)
      }
  }

  def parseFormatKeys[_: P]: P[Seq[VcfFormatHeader]] = formatKey.rep(1, sep=":")
}

sealed trait GenotypeParser {
  def header: VcfHeader

  def gtValue[_: P]: P[SI] = (Index ~ ("0".! | (CharIn("1-0") ~ CharIn("0-9").rep).!) ~ Index).si
  def gt[_: P](alleles: AlleleSet): P[Allele] = ".".!.map(_ => NoCallAllele) | gtValue.flatMap { v =>
    val index = v.value.toInt
    alleles.get(index) match {
      case None         => v.fail(s"genotype index `$index` is out of range: [0, ${alleles.size}]")
      case Some(allele) => success(allele)
    }
  }
  def parseGT[_: P](alleles: AlleleSet): P[Seq[Allele]] = {
    (gt(alleles).rep(min=1, sep="/") | gt(alleles).rep(min=1, sep="|"))
  }

//  def sampleValue[_: P]: P[Any] = (Index ~ (".".map | CharPred(_ != ':').rep(1)).! ~ Index).si.flatMap { key =>
//    if (key.value)
//
//  }
  def sampleValues[_: P]: P[Seq[String]] = (Start ~ sampleValue.rep(0, sep=":") ~ End)


  def parseGenotype[_: P](formatKeys: Seq[VcfFormatHeader], alleles: AlleleSet): P[Genotype] = {

    formatKeys match {
      case Seq(gtKey, tail @ _*) if gtKey.id == "GT" =>
        parseGT(alleles)

      case _                                         =>
    }
  }

}

class VariantParser(val header: VcfHeader)
  extends SequenceMetadataParser
    with PositionParser
    with IdentifierParser
    with RefAlleleParser
    with AltAlleleParser
    with QualParser
    with FilterParser
    with InfoParser
    with FormatKeysParser {

  private val numSamples: Int = header.samples.length

  private def parseChromPosAndIdentifier[_: P]: P[(SequenceMetadata, Int, Seq[SI])] = {
    Start ~/
      parseSequenceMetadata ~/ // CHROM
      "\t" ~/ parsePosition ~/ // POS
      "\t" ~/ parseIdentifier ~/  // ID
      "\t"
  }

  def variant[_: P]: P[Variant] = {
    parseChromPosAndIdentifier.flatMap { case (metadata, pos, id) =>
      parseRefAllele(metadata=metadata, pos=pos).map { ref => (metadata, pos, id, ref) } ~/ "\t"
    }.flatMap { case (metadata, pos, id, ref) =>
      (
        parseAltAlleles(ref=ref) ~/ // ALT
        "\t" ~/ parseQual ~/        // QUAL
        "\t" ~/ parseFilter ~/      // FILTER
        "\t"
      ).map { case (alleleSet, qual, filters) =>
        (metadata, pos, id, ref, alleleSet, qual, filters)
      }
    }.flatMap { case (metadata, pos, id, ref, alleleSet, qual, filters) =>
      // INFO
      parseInfos(alleleSet).map { infos =>
        (metadata, pos, id, ref, alleleSet, qual, filters, infos)
      }
    }.flatMap { case (metadata, pos, id, ref, alleleSet, qual, filters, infos) =>
      if (this.numSamples == 0) {
        "" ~ End.map { _ =>
          (metadata, pos, id, ref, alleleSet, qual, filters, infos, Seq.empty[Genotype])
        }
      }
      else {
        ("\t" ~/ parseFormatKeys).flatMap { formats =>
          // FIXME
          success(metadata, pos, id, ref, alleleSet, qual, filters, infos, Seq.empty[Genotype])
          // TODO: parse a sample!
        }
      }
    }.map { case (metadata, pos, id, ref, alleleSet, qual, filters, infos, genotypes) =>
      new Variant(
        chrom     = metadata.name,
        pos       = pos,
        id        = id.map(_.value),
        alleles   = alleleSet,
        qual      = qual,
        filters   = filters.map(_.id).toSet,
        attrs     = infos,
        genotypes = genotypes.map(genotype => genotype.sample -> genotype).toMap
      )
    }
  }
}
