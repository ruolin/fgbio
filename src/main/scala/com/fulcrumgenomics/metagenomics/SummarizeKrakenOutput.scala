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

import java.text.DecimalFormat

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.metagenomics.NcbiTaxonomy.TaxonId
import com.fulcrumgenomics.util.{Io, ProgressLogger, SimpleCounter}
import dagr.commons.util.LazyLogging
import dagr.sopt._

import scala.annotation.tailrec
import scala.collection.mutable.ListBuffer

@clp(group=ClpGroups.Metagenomics, description=
"""
  |Summarizes the output of running Kraken[1] on a dataset. Summarization takes the
  |form of pushing low-level taxa assignments "up" the tree until the assignment is
  |at or above a given rank (e.g. genus), and pushing high-level taxa assignments
  |down until at the given rank.
  |
  |[1] Kraken: http://ccb.jhu.edu/software/kraken/
""")
class SummarizeKrakenOutput
( @arg(flag="i", doc="The kraken output file.")     val input: FilePath,
  @arg(flag="o", doc="Output summary report file.") val output: FilePath,
  @arg(flag="d", doc="The kraken database used.")   val database: DirPath,
  @arg(flag="r", doc="The taxanomic rank to summarize at.") val rank: Rank = Rank.genus,
  @arg(flag="m", doc="Minimum number of hits to support rank assignment.") val minimumHits: Int = 2,
  @arg(flag="p", doc="Push rank assignments 'up' until all kmer hits are in-tree.")  val push: Boolean = true,
  @arg(flag="f", doc="Filter kmer hits to non-orgamism parts of the taxonomy tree.") val filterOther: Boolean = true
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertListable(database)
  Io.assertCanWriteFile(output)

  private val taxaDir = database.resolve("taxonomy")
  if (!taxaDir.toFile.exists()) fail("Database does not have 'taxonomy' directory.")
  logger.info("Loading taxonomy information.")
  private val taxa = NcbiTaxonomy(names=taxaDir.resolve("names.dmp"), nodes=taxaDir.resolve("nodes.dmp"))

  // Some useful parts of the taxonomy tree
  private val RootTaxonId = 1
  private val UnassignedTaxonId = 0
  private val OtherSequencesId = 28384 // taxon id for "other" sequences like artificial sequences and non-organisms
  private val OtherSequencesTaxon = this.taxa(OtherSequencesId)

  override def execute(): Unit = {
    logger.info("Loading kraken results.")
    val progress = new ProgressLogger(logger, noun="records")
    val counter = new SimpleCounter[TaxonId]
    new KrakenReader(input).foreach { rec =>
      if (this.taxa.get(rec.taxonId).isDefined) {
        val filtered = if (filterOther) filterOther(rec) else rec
        val pushed = if (push) push(filtered) else rec
        val supported = hitSupportedRank(pushed)
        counter.count(supported.taxonId)
      }

      progress.record()
    }

    logger.info(s"Raw     = Queries: ${counter.total}; Hits: ${counter.total-counter.countOf(0)}; Taxa: ${counter.size-1}.")

    // Roll things down from higher levels to lower levels
    val tree = buildSummaryTree(counter, taxa)
    tree.rolldown(rank)
    tree.rollup(rank)

    logger.info(s"Summary = Queries: ${counter.total}; Hits: ${tree.total}; Taxa: ${tree.count(_.count > 0)}.")

    writeReport(tree, counter.countOf(0))
  }

  /**
    * Filters out hits to 'other' taxonomy entries like artificial sequences etc.
    * This will also have the side effect of filtering out all hits to retired
    * taxonomy IDs that are present in sequences, but not in the taxonomy db anymore.
    */
  def filterOther(rec: KrakenResult): KrakenResult = {
    val filteredCounts = rec.kmerCounts.filter { case (id, n) =>
        id == UnassignedTaxonId || taxa.get(id).exists(t => !t.descendsFrom(this.OtherSequencesTaxon))
    }

    rec.copy(kmerCounts = filteredCounts)
  }

  @tailrec
  private[metagenomics] final def push(rec: KrakenResult): KrakenResult = {
    if (rec.taxonId == UnassignedTaxonId || rec.taxonId == RootTaxonId) {
      rec
    }
    else {
      val taxon = taxa(rec.taxonId)
      val ok = rec.kmerCounts.forall {
        case (UnassignedTaxonId, n) => true
        case (RootTaxonId, n)       => true
        case (id, n) =>
          val other = taxa(id)
          other.id == taxon.id || other.descendsFrom(taxon) || taxon.descendsFrom(other)
      }

      (ok, taxon.parent) match {
        case (true, _)        => rec.copy(taxonId = taxon.id)
        case (false, None)    => rec.copy(taxonId = UnassignedTaxonId)
        case (false, Some(p)) => push(rec.copy(taxonId=p.id))
      }
    }
  }


  /**
    * Looks at the kmer assignments for a kraken result and pushes the assignment up the tree
    * until at least this.hits kmers have matched.
    * @param rec the kraken result
    * @return the rank that the result should be counted at
    */
  private[metagenomics] def hitSupportedRank(rec: KrakenResult): KrakenResult = {
    if (rec.taxonId == UnassignedTaxonId || rec.kmerCounts.getOrElse(rec.taxonId, 0) >= this.minimumHits) {
      rec
    }
    else {
      val counts = rec.kmerCounts
      var taxon   = taxa(rec.taxonId)
      var support = subtreeCount(rec, taxa)
      while (support < this.minimumHits && taxon.parent.isDefined) {
        taxon = taxon.parent.getOrElse(unreachable("taxon.parent.isDefined checked directly above."))
        support += counts.getOrElse(taxon.id, 0)
      }

      if (support >= this.minimumHits) rec.copy(taxonId=taxon.id)
      else rec.copy(taxonId=UnassignedTaxonId)
    }
  }

  /** Gets the count of kmers that matched either the assigned taxon, or descendants in the taxon tree. */
  private[metagenomics] def subtreeCount(rec: KrakenResult, taxa: NcbiTaxonomy): Int = rec.taxonId match {
    case UnassignedTaxonId => 0
    case taxonId =>
      var count = 0
      val taxon = taxa(taxonId)
      rec.kmerCounts.foreach { case(id, n) => if (id != UnassignedTaxonId && taxa.get(id).forall(_.descendsFrom(taxon))) count += n }
      count
  }

  /** A little class that is similar to Taxon, but carries around a count. */
  private class CountedTaxon(val id: Int, val name: String, val parent: Option[CountedTaxon], val rank: Option[Rank], var count: Double=0)
    extends Iterable[CountedTaxon] {
    val children = ListBuffer[CountedTaxon]()

    /** Gets the total count of this node and all it's descendant nodes. */
    def total: Double = count + children.map(_.total).sum

    /** Eliminates children below a certain rank. */
    def rollup(lowestRank: Rank) : Unit = {
      this.children.foreach { child =>
        child.rollup(lowestRank)
        if (child.rank.exists(lowestRank.isAbove)) {
          this.count += child.count
          child.count = 0
        }
      }
    }

    def rolldown(rank: Rank) : Unit = {
      if (this.rank.forall(_.isAbove(rank))) {
        if (this.count > 0) {
          val kidCounts = this.children.map(_.total)
          val kidTotal  = kidCounts.sum
          val kidsWithTotals = this.children.zip(kidCounts).filter(_._2 > 0)

          if (kidTotal > 0 && kidsWithTotals.map(_._1).forall(kid => kid.rank.forall(r => r.isAtOrAbove(rank)))) {
            kidsWithTotals.foreach { case (kid, n) =>
              val proportion = n / kidTotal
              kid.count += this.count * proportion
            }

            this.count = 0
          }
        }

        this.children.foreach(_.rolldown(rank))
      }
    }

    /** Provides an iterator over all CountedTaxons starting from this node. */
    override def iterator: Iterator[CountedTaxon] = Iterator(this) ++ this.children.iterator.flatMap(_.iterator)
  }

  /** Converts the tree of Taxons into a smaller tree of CountedTaxons containing only those branches that
    * lead to nodes with counts.
    */
  private def buildSummaryTree(counter: SimpleCounter[TaxonId], taxa: NcbiTaxonomy) : CountedTaxon = {
    // The set of tax ids that have counts, or have children with counts
    val ids = counter.flatMap(kv => taxa.get(kv._1)).flatMap(t => t :: t.ancestors.toList).map(_.id).toSet

    /** Recursive function to build out children nodes of a given node. */
    def buildChildren(parent: Taxon, countedParent: CountedTaxon) : Unit = {
      parent.children.foreach { child =>
        if (ids.contains(child.id)) {
          val counted = new CountedTaxon(id=child.id, name=child.name, parent=Some(countedParent), rank=child.rank, count=counter.countOf(child.id))
          buildChildren(child, counted)
          countedParent.children += counted
        }
      }
    }

    val root    = taxa.root
    val counted = new CountedTaxon(id=root.id, name=root.name, parent=None, rank=root.rank, count=counter.countOf(root.id))
    buildChildren(root, counted)
    counted
  }

  def writeReport(root: CountedTaxon, unassigned: Double): Unit = {
    case class ReportLine(id: Int, name: String, rank: String, count: Double)
    val ffmt = new DecimalFormat("0.000000")
    val nfmt = new DecimalFormat("0.00")
    val totalAssigned = root.total
    val total = totalAssigned + unassigned
    val lines = ReportLine(id=0, name="unassigned", rank="n/a", count=unassigned) ::
      root.filter(_.count > 0).map(t => ReportLine(id=t.id, name=t.name, rank=t.rank.getOrElse("n/a").toString, count=t.count)).toList
    val out = Io.toWriter(output)
    out.write(Seq("taxon_id", "name", "rank", "hits", "fraction_of_total", "fraction_of_assigned").mkString("\t"))
    out.newLine()

    // TODO: cumulative percentage?

    lines.sortBy(- _.count).foreach { line =>
      val fields = Seq(line.id, line.name, line.rank,
                       nfmt.format(line.count),
                       ffmt.format(line.count / total),
                       if (line.id == 0) 0 else  ffmt.format(line.count / totalAssigned)
      )
      out.write(fields.mkString("\t"))
      out.newLine()
    }

    out.close()
  }
}
