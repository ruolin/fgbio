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

package com.fulcrumgenomics.bam.api


import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Cigar
import htsjdk.samtools
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools._
import htsjdk.samtools.util.CoordMath

// TODO long-term: methods for 5'end, 3' end, unclipped 5' end and unclipped 3' end
// TODO long-term: replacement for alignment blocks?

/**
  * Trait that extends [[SAMRecord]] so that it can override methods and call super.  Is separate from
  * [[SamRecord]] so that [[SamRecord]] does not export the [[SAMRecord]] API.
  *
  * Currently this class is only necessary because a) we define our own [[Cigar]] class and cache the
  * parsed value in [[SamRecord]] and b) we use HTSJDK utilities that may set the cigar using the
  * [[SAMRecord]] API methods `setCigar()` and `setCigarString()`.
  */
private[api] trait SamRecordIntermediate extends SAMRecord {
  final override def setCigar(cigar: samtools.Cigar): Unit = {
    val s = if (cigar == null) null else TextCigarCodec.encode(cigar)
    setCigarString(s)
  }

  override def setCigarString(cigar: String): Unit = { super.setCigarString(cigar);  cigarChanged(cigar) }
  protected[api] def setCigarStringNoNotify(cigar: String): Unit = super.setCigarString(cigar)
  protected[api] def cigarChanged(cigar: String): Unit = {}
}

/** Class that is used to provide a nice API to transient attributes in the SamRecord. */
class TransientAttrs(private val rec: SamRecord) {
  def apply[A](key: Any): A = rec.asSam.getTransientAttribute(key).asInstanceOf[A]
  def update(key: Any, value: Any): Unit = {
    if (value == null) rec.asSam.removeTransientAttribute(key) else rec.asSam.setTransientAttribute(key, value)
  }
  def get[A](key: Any): Option[A] = Option(apply(key))
  def getOrElse[A](key: Any, default: => A): A = rec.asSam.getTransientAttribute(key) match {
    case null => default
    case value => value.asInstanceOf[A]
  }
  def getOrElseUpdate[A](key: Any, default: => A): A = rec.asSam.getTransientAttribute(key) match {
    case null  =>
      val value: A = default
      update(key, value)
      value
    case value => value.asInstanceOf[A]
  }
}

/**
  * A trait that fgbio uses as a replacement for [[SAMRecord]].  The trait is self-typed as a
  * [[SamRecordIntermediate]] which is a sub-class of SAMRecord.  It is done this wasy so that
  * a) we can access superclass methods via [[SamRecordIntermediate]] but that self-typing here
  * instead of extending hides the [[SAMRecord]] API from users of the class.  The result is
  * always a [[SAMRecord]] but isn't seen as such without casting.
  */
trait SamRecord {
  this: SamRecordIntermediate =>

  private var _cigar: Cigar = _

  ////////////////////////////////////////////////////////////////////////////
  // Wrapper methods.
  ////////////////////////////////////////////////////////////////////////////

  @inline final def name: String = getReadName
  @inline final def name_=(name: String): Unit = setReadName(name)

  @inline final def flags: Int = getFlags

  @inline final def paired: Boolean = getReadPairedFlag
  @inline final def paired_=(paired: Boolean):Unit = setReadPairedFlag(paired)

  @inline final def unpaired: Boolean = !paired
  @inline final def unpaired_=(unpaired: Boolean): Unit = this.paired = !unpaired

  @inline final def properlyPaired: Boolean = getProperPairFlag
  @inline final def properlyPaired_=(paired: Boolean):Unit = setProperPairFlag(paired)

  @inline final def mapped: Boolean = !getReadUnmappedFlag
  @inline final def mapped_=(mapped: Boolean):Unit = setReadUnmappedFlag(!mapped)
  @inline final def unmapped: Boolean = getReadUnmappedFlag
  @inline final def unmapped_=(unmapped: Boolean):Unit = setReadUnmappedFlag(unmapped)

  @inline final def mateMapped: Boolean = !getMateUnmappedFlag
  @inline final def mateMapped_=(mapped: Boolean):Unit = setMateUnmappedFlag(!mapped)
  @inline final def mateUnmapped: Boolean = getMateUnmappedFlag
  @inline final def mateUnmapped_=(unmapped: Boolean):Unit = setMateUnmappedFlag(unmapped)

  @inline final def positiveStrand: Boolean = !getReadNegativeStrandFlag
  @inline final def positiveStrand_=(positive: Boolean):Unit = setReadNegativeStrandFlag(!positive)
  @inline final def negativeStrand: Boolean = getReadNegativeStrandFlag
  @inline final def negativeStrand_=(negative: Boolean):Unit = setReadNegativeStrandFlag(negative)

  @inline final def matePositiveStrand: Boolean = !getMateNegativeStrandFlag
  @inline final def matePositiveStrand_=(positive: Boolean):Unit = setMateNegativeStrandFlag(!positive)
  @inline final def mateNegativeStrand: Boolean = getMateNegativeStrandFlag
  @inline final def mateNegativeStrand_=(negative: Boolean):Unit = setMateNegativeStrandFlag(negative)

  @inline final def firstOfPair: Boolean = getFirstOfPairFlag
  @inline final def firstOfPair_=(first: Boolean):Unit = setFirstOfPairFlag(first)

  @inline final def secondOfPair: Boolean = getSecondOfPairFlag
  @inline final def secondOfPair_=(second: Boolean):Unit = setSecondOfPairFlag(second)

  @inline final def secondary: Boolean = isSecondaryAlignment
  @inline final def secondary_=(secondary: Boolean):Unit = setSecondaryAlignment(secondary)
  @inline final def primary: Boolean = !this.secondary
  @inline final def primary_=(primary: Boolean):Unit = this.secondary = !primary

  @inline final def pf: Boolean = !getReadFailsVendorQualityCheckFlag
  @inline final def pf_=(pf: Boolean):Unit = setReadFailsVendorQualityCheckFlag(!pf)

  @inline final def duplicate: Boolean = getDuplicateReadFlag
  @inline final def duplicate_=(dupe: Boolean):Unit = setDuplicateReadFlag(dupe)

  @inline final def supplementary: Boolean = getSupplementaryAlignmentFlag
  @inline final def supplementary_=(supplementary: Boolean):Unit = setSupplementaryAlignmentFlag(supplementary)

  @inline final def refName: String = getReferenceName
  @inline final def refName_=(name: String):Unit = setReferenceName(name)

  @inline final def refIndex: Int = getReferenceIndex
  @inline final def refIndex_=(index: Int):Unit = setReferenceIndex(index)

  @inline final def start: Int = getAlignmentStart
  @inline final def start_=(s: Int):Unit = setAlignmentStart(s)

  @inline final def end: Int = if (unmapped) SAMRecord.NO_ALIGNMENT_START else start + cigar.lengthOnTarget - 1

  @inline final def unclippedStart: Int = if (unmapped) SAMRecord.NO_ALIGNMENT_START else getUnclippedStart
  @inline final def unclippedEnd  : Int = if (unmapped) SAMRecord.NO_ALIGNMENT_START else getUnclippedEnd

  @inline final def mapq: Int = getMappingQuality
  @inline final def mapq_=(q: Int):Unit = setMappingQuality(q)

  // The cigar handling code here is rather ugly because we can't _both_ override setCigar/setCigarString
  // and call super.setCigar/setCigarString in a self-typed trait :(  So we need to leave those methods
  // alone, but also detect when they've been called and reset our cigar representation!
  @inline final def cigar: Cigar = { if (_cigar == null) _cigar = Cigar.fromSam(getCigarString);  _cigar }
  @inline final def cigar_=(cig: String): Unit = { this._cigar = null; setCigarStringNoNotify(cig) }
  @inline final def cigar_=(cig: Cigar): Unit  = { this._cigar = cig; setCigarStringNoNotify(if (cig.isEmpty) null else cig.toString()) }
  @inline final override def cigarChanged(cigar: String): Unit = { this._cigar = null } // null out the cigar; lazily reconstruct if asked for again

  @inline final def mateRefName: String = getMateReferenceName
  @inline final def mateRefName_=(name: String):Unit = setMateReferenceName(name)

  @inline final def mateRefIndex: Int = getMateReferenceIndex
  @inline final def mateRefIndex_=(index: Int):Unit = setMateReferenceIndex(index)

  @inline final def mateStart: Int = getMateAlignmentStart
  @inline final def mateStart_=(s: Int):Unit = setMateAlignmentStart(s)

  @inline final def mateCigar: Option[Cigar] = {
    require(paired && mateMapped, "Cannot get a mate cigar on read without a mapped mate.")
    get[String]("MC").map(Cigar.apply)
  }
  @inline final def mateEnd: Option[Int] = {
    require(paired && mateMapped, "Cannot get mate end position on read without a mapped mate.")
    get[String]("MC").map(cig => mateStart + Cigar(cig).lengthOnTarget - 1)
  }
  @inline final def mateUnclippedStart: Option[Int] = {
    require(paired && mateMapped, "Cannot get mate unclipped start position on read without a mapped mate.")
    get[String]("MC").map(cig => mateStart - Cigar(cig).iterator.takeWhile(_.operator.isClipping).map(_.length).sum)
  }
  @inline final def mateUnclippedEnd: Option[Int] = {
    require(paired && mateMapped, "Cannot get mate unclipped end position on read without a mapped mate.")
    get[String]("MC").map { cig =>
      val cigar = Cigar(cig)
      mateStart + cigar.lengthOnTarget - 1 + cigar.reverseIterator.takeWhile(_.operator.isClipping).map(_.length).sum
    }
  }
  @inline final def matesOverlap: Option[Boolean] = {
    require(mapped && paired && mateMapped, "Cannot determine if mates overlap without paired mates that are both mapped.")
    if (refIndex != mateRefIndex) Some(false)
    else if (mateStart > end) Some(false)
    else if (mateStart >= start && mateStart <= end) Some(true)
    else mateEnd.map(mateEnd => CoordMath.overlaps(start, end, mateStart, mateEnd))
  }

  @inline final def insertSize: Int = getInferredInsertSize
  @inline final def insertSize_=(s: Int):Unit = setInferredInsertSize(s)

  @inline final def bases: Array[Byte] = getReadBases
  @inline final def basesString: String = getReadString
  @inline final def bases_=(bs: Array[Byte]): Unit = setReadBases(bs)
  @inline final def bases_=(bs: String): Unit = setReadString(bs)

  @inline final def quals: Array[Byte] = getBaseQualities
  @inline final def qualsString: String = getBaseQualityString
  @inline final def quals_=(qs: Array[Byte]): Unit = setBaseQualities(qs)
  @inline final def quals_=(qs: String): Unit = setBaseQualityString(qs)

  @inline final def length: Int = getReadLength

  @inline final def header: SAMFileHeader = getHeader
  @inline final def header_=(header: SAMFileHeader) = setHeader(header)
  @inline final def readGroup: SAMReadGroupRecord = getReadGroup

  // Use apply/update for tag attributes
  @inline final def apply[A](name: String): A                    = getAttribute(name).asInstanceOf[A]
  @inline final def get[A](name: String): Option[A]              = Option(apply(name))
  @inline final def update(name: String, value: Any): Unit       = setAttribute(name, value)
  @inline final def attributes: Map[String,Any]                  = getAttributes.map(x => x.tag -> x.value).toMap
  @inline final def getOrElse[A](name: String, default: => A): A = get(name).getOrElse(default)
  @inline final def contains(name: String): Boolean              = hasAttribute(name)
  @inline final def remove(name: String): Unit                   = setAttribute(name, null)

  // transient attributes
  @inline final val transientAttrs: TransientAttrs = new TransientAttrs(this)

  /** Returns the 1-based reference position for the 1-based position in the read, or [[None]] if the reference does
    * not map at the read position (e.g. insertion).
    *
    * @param pos the 1-based read position to query
    * @param returnLastBaseIfInserted if the reference is an insertion, true to returned the previous reference base,
    *                                 false to return None
    * */
  @inline final def refPosAtReadPos(pos: Int, returnLastBaseIfInserted: Boolean = false): Option[Int] = if (pos <= 0) None else {
    var readPos      = 1
    var refPos       = this.start
    var elementIndex = 0
    val elems        = this.cigar.elems

    def continue(): Boolean = if (elementIndex >= elems.length) false else {
      val elem = elems(elementIndex)
      if (pos > CoordMath.getEnd(readPos, elem.lengthOnQuery)) true // current cigar element is before the desired position
      else {
        if (pos < readPos) readPos = 0 // past the point in the read
        else if (elem.operator == CigarOperator.INSERTION) {
          if (!returnLastBaseIfInserted) refPos = 0
          else refPos -= 1
        } // in an insertion, no reference position
        else refPos += pos - readPos  // get the offset
        false
      }
    }

    while (continue()) {
      val elem      = elems(elementIndex)
      refPos       += elem.lengthOnTarget
      readPos      += elem.lengthOnQuery
      elementIndex += 1
    }

    if (this.start <= refPos && readPos <= this.end) Some(refPos) else None
  }


  /** Returns the 1-based position into the read's bases where the position is mapped either as a match or mismatch,
    * [[None]] if the read does not map the position.
    *
    * @param pos the 1-based reference position to query
    * @param returnLastBaseIfDeleted if the reference is a deletion, true to returned the previous read base, false to
    *                                return None
    * */
  @inline final def readPosAtRefPos(pos: Int, returnLastBaseIfDeleted: Boolean = false): Option[Int] = {
    if (this.unmapped || pos < this.start || this.end < pos) None
    else { // overlaps
      var readPos      = 0
      var refPos       = this.start
      var elementIndex = 0
      val elems        = this.cigar.elems

      def continue(): Boolean = if (elementIndex >= elems.length) false else {
        val elem = elems(elementIndex)
        if (pos > CoordMath.getEnd(refPos, elem.lengthOnTarget)) true // current cigar element is before the desired position
        else { // overlaps!
          if (elem.operator == CigarOperator.DELETION) { // deletion
            if (!returnLastBaseIfDeleted) readPos = 0 // don't return a read position, so zero out the current read position
          } else {
            readPos += pos - refPos + 1 // get the offset
          }
          false
        }
      }

      while (continue()) {
        val elem      = elems(elementIndex)
        refPos       += elem.lengthOnTarget
        readPos      += elem.lengthOnQuery
        elementIndex += 1
      }

      if (0 < readPos && readPos <= this.length) Some(readPos) else None
    }
  }


  ////////////////////////////////////////////////////////////////////////////
  // Non-wrapper methods.
  ////////////////////////////////////////////////////////////////////////////

  /** Returns a string that is useful to identify a SamRecord, mostly for testing and error messages. */
  def id: String = {
    val builder = new StringBuilder
    builder.append(name)
    builder.append(if (paired && secondOfPair) "/2" else "/1")
    if (secondary) builder.append(":sec")
    if (supplementary) builder.append(":sup")
    builder.toString()
  }

  /** Returns this record as a SAMRecord. */
  def asSam: SAMRecord = this.asInstanceOf[SAMRecord]

  /** Gets the PairOrientation of the record. */
  def pairOrientation: PairOrientation = SamPairUtil.getPairOrientation(this.asSam)

  /** Returns true if the read is mapped in an FR pair, false otherwise. */
  def isFrPair: Boolean = {
    paired &&
      mapped &&
      mateMapped &&
      refIndex == mateRefIndex &&
      SamPairUtil.getPairOrientation(this) == PairOrientation.FR
  }

  /** Clone method that does a "reasonably deep" clone. The bases and quals are cloned as is the attributes map,
    * though not the values in the attributes map. */
  override def clone(): SamRecord = {
    val r = super.clone().asInstanceOf[SamRecord]
    r.bases = this.bases.clone()
    r.quals = this.quals.clone()
    r
  }
}

object SamRecord {

  /** The mapping quality when not available. */
  val UnknownMappingQuality: Int = SAMRecord.UNKNOWN_MAPPING_QUALITY

  /** The mapping quality for an unmapped or multi-mapped read */
  val ZeroMappingQuality: Int = SAMRecord.NO_MAPPING_QUALITY

  /** The reference name for an unmapped read. */
  val UnmappedReferenceName: String = SAMRecord.NO_ALIGNMENT_REFERENCE_NAME

  /** The reference index for an unmapped read. */
  val UnmappedReferenceIndex: Int = SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX

  /** The cigar string for an unmapped read. */
  val UnmappedCigarString: String = SAMRecord.NO_ALIGNMENT_CIGAR

  /** The reference position for an unmapped read. */
  val UnmappedStart: Int = SAMRecord.NO_ALIGNMENT_START

  /** The value for the base string when no bases are present. */
  val MissingBases: String = SAMRecord.NULL_SEQUENCE_STRING

  /** The value for the base quality string when no base qualities are present. */
  val MissingQuals: String = SAMRecord.NULL_QUALS_STRING

  /** The maximum insert size that can be stored in a [[SamRecord]]. */
  val MaximumInsertSize: Int = SAMRecord.MAX_INSERT_SIZE

  /** SAMRecord with mixin to add behaviour. */
  private final class EnhancedSamRecord(header: SAMFileHeader) extends SAMRecord(header) with SamRecordIntermediate with SamRecord

  /** BAMRecord with mixin to add behaviour. */
  private final class EnhancedBamRecord(header: SAMFileHeader,
                                        refIndex: Int,
                                        alignmentStart: Int,
                                        readNameLength: Short,
                                        mapq: Short,
                                        indexingBin: Int,
                                        cigarLen: Int,
                                        flags: Int,
                                        readLen: Int,
                                        mateRefIndex: Int,
                                        mateAlignmentStart: Int,
                                        insertSize: Int,
                                        variableLengthBlock: Array[Byte]) extends
    BAMRecord(header, refIndex, alignmentStart, readNameLength, mapq, indexingBin, cigarLen,
      flags, readLen, mateRefIndex, mateAlignmentStart, insertSize, variableLengthBlock) with SamRecordIntermediate with SamRecord

  /** Factory method for new [[SamRecord]] instances. */
  def apply(header: SAMFileHeader): SamRecord = new EnhancedSamRecord(header)

  /**
    * Singleton implementation of SAMRecordFactory to return instances of [[htsjdk.samtools.SAMRecord]] and
    * [[htsjdk.samtools.BAMRecord]] that have been enhanced by mixing in [[SamRecord]].
    */
  object Factory extends SAMRecordFactory {
    /** Factory implementation method to create SAMRecords. */
    override def createSAMRecord(header: SAMFileHeader): SAMRecord = new EnhancedSamRecord(header)

    /** Factory implementation method to create BAMRecords. */
    override def createBAMRecord(header: SAMFileHeader,
                                 refIndex: Int,
                                 alignmentStart: Int,
                                 readNameLength: Short,
                                 mapq: Short,
                                 indexingBin: Int,
                                 cigarLen: Int,
                                 flags: Int,
                                 readLen: Int,
                                 mateRefIndex: Int,
                                 mateAlignmentStart: Int,
                                 insertSize: Int,
                                 variableLengthBlock: Array[Byte]): BAMRecord =
      new EnhancedBamRecord(header, refIndex, alignmentStart, readNameLength, mapq, indexingBin, cigarLen,
        flags, readLen, mateRefIndex, mateAlignmentStart, insertSize, variableLengthBlock)
  }
}
