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

import java.nio.file.Path

import com.fulcrumgenomics.metagenomics.NcbiTaxonomy.TaxonId
import com.fulcrumgenomics.util.Io

import scala.annotation.tailrec
import scala.collection.mutable

/**
  * Companion object for NcbiTaxonomy that provides methods for parsing NCBI taxonomy
  * information, and a handful of useful types/constants.
  */
object NcbiTaxonomy {
  type TaxonId = Int
  private val SplitPattern = "\t\\|\t?".r
  private val ScientificNameType = "scientific name"

  /** Loads up the NCBI taxonomy from the names.dmp and nodes.dmp files. */
  def apply(names: Path, nodes: Path): NcbiTaxonomy = {
    val idToTaxon = mutable.HashMap[TaxonId,Taxon]()
    open(names).filter(fs => fs(3) == ScientificNameType).foreach { fs =>
      val id = fs(0).toInt
      val name = fs(1).trim
      val taxon = new Taxon(id=id, name=name)
      idToTaxon(id) = taxon
    }

    open(nodes).foreach { fs =>
      val id       = fs(0).toInt
      val parentId = fs(1).toInt
      val rank     = if (fs(2) == "no rank") None else Some(Rank.of(fs(2)))

      val taxon  = idToTaxon(id)
      val parent = idToTaxon(parentId)

      taxon.rank = rank
      if (taxon != parent) {
        taxon.parent = Some(parent)
        parent._children.add(taxon)
      }
    }

    new NcbiTaxonomy(idToTaxon.toMap)
  }

  /** Opens an NCBI Tax dump file that uses tab-pipe-tab as a separator and returns it as an iterator lines as fields. */
  private def open(dump: Path) : Iterator[Array[String]] = Io.toSource(dump).getLines().map(SplitPattern.split)
}


/**
  * A class representing the entire NCBI taxonomy, or any portion of it loaded.
  */
class NcbiTaxonomy(private val taxaById: Map[TaxonId,Taxon]) {
  /** Returns the root entry if present in the taxonomy. */
  val root: Taxon = taxaById(1)

  /** Looks up a Taxon by ID. */
  def apply(id: TaxonId): Taxon = this.taxaById(id)

  /** Looks up a Taxon by ID, returning an option to handle possibly absent values. */
  def get(id: TaxonId): Option[Taxon] = this.taxaById.get(id)

  /** Provides an iterator over the set of TaxonIds stored in the taxonomy. */
  def ids: Iterable[TaxonId] = this.taxaById.keys

  /** Provides an iterator over the set of Taxa stored in the taxonomy. */
  def taxa: Iterable[Taxon] = this.taxaById.values
}

/**
  * An individual Taxon within the Taxonomy.
  * @param id the ID of the taxon
  * @param name the unique scientific name assigned to the taxon
  * @param rank the rank within the taxonomy if one is assigned
  * @param parent the parent taxon, if there is one
  */
class Taxon(val id: TaxonId, val name: String, var rank: Option[Rank]=None, var parent: Option[Taxon]=None) {
  private[metagenomics] val _children: mutable.Set[Taxon] = mutable.Set()

  /** Finds the first ancestor that is at or above the provided rank. Since not all branches of
    * the taxonomy descend through all ranks, it is possible to return a parent at a higher rank
    * than requested.
    *
    * @param lowestRank the lowest desired taxanomic rank of the parent
    * @return Some(Taxon) if there exists a parent with a defined rank that is at least as high as
    *         lowestRank, false otherwise.
    */
  def findAncestor(lowestRank: Rank): Option[Taxon] = (this.rank, this.parent) match {
    case (Some(r), _) if r.isAtOrAbove(lowestRank) => Some(this)
    case (_, Some(p))                              => p.findAncestor(lowestRank)
    case _                                         => None
  }

  /** Returns an iterator over all ancestors of this node. */
  def ancestors: Iterator[Taxon] = parent match {
    case None    => Iterator.empty
    case Some(p) => Iterator(p) ++ p.ancestors
  }

  /** Returns true if the other taxon is an ancestor of this taxon and false otherwise. */
  def descendsFrom(other: Taxon): Boolean = {
    var taxon = this.parent
    while (taxon.isDefined && taxon.get.id != other.id) taxon = taxon.get.parent
    taxon.isDefined
  }

  /** Returns an iterator over all direct children of this taxon. */
  def children: Iterator[Taxon] = this._children.iterator

  /** Returns all descendants of this node all the way down to the leaves. */
  def descendants: Iterator[Taxon] = this.children ++ this.children.flatMap(_.descendants)

  /** Returns this node and all of it's descendants. */
  def subtree: Iterator[Taxon] = Iterator(this) ++ descendants

  /** Returns the 'level' in the tree, where the root node is level 0, children of root are level 1 and so on. */
  def level: Int = parent match {
    case None    => 0
    case Some(p) => 1 + p.level
  }

  /** Equality based solely on the ID. */
  override def equals(other: Any): Boolean = other match {
    case that: Taxon => id == that.id
    case _           => false
  }

  /** Hashcode based on ID only. */
  override def hashCode(): Int = id.hashCode()

  override def toString: String = render()

  @tailrec
  private def render(taxon: Option[Taxon]=Some(this), sofar: StringBuilder=new StringBuilder): String = taxon match {
    case None    => sofar.toString()
    case Some(t) =>
      if (sofar.nonEmpty) sofar.append(" -> ")
      sofar.append("Taxon(id=").append(t.id).append(", name=").append(t.name).append(")")
      render(taxon=t.parent, sofar)
  }
}

