/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

import java.io.PrintStream
import java.nio.file.Paths
import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Aligner
import com.fulcrumgenomics.alignment.Mode.Global
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.fastq.{FastqSource, QualityEncoding}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric, ProgressLogger, Rscript}
import htsjdk.samtools.CigarOperator

import math.{max, min}
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer

object CheckDropSeqData {
  private val PlotScript = "com/fulcrumgenomics/personal/tfenne/checkDropSeqData.R"

  /** Summary metrics. */
  case class DropSeqSummaryMetrics(sample: String,
                                   expected_cells: Long,
                                   data_type: String,
                                   read_count: Long,
                                   cbc_count: Long,
                                   umi_count: Long,
                                   median_umis_top_500_cells: Double,
                                   mean_umis_top_500_cells: Double,
                                   median_umis_top_1000_cells: Double,
                                   mean_umis_top_1000_cells: Double,
                                   median_umis_top_1500_cells: Double,
                                   mean_umis_top_1500_cells: Double,
                                   median_umis_top_5000_cells: Double,
                                   mean_umis_top_5000_cells: Double,
                                   median_umis_top_10000_cells: Double,
                                   mean_umis_top_10000_cells: Double,
                                   fraction_cbcs_with_250_umis: Double,
                                   fraction_cbcs_with_500_umis: Double,
                                   fraction_cbcs_with_1000_umis: Double
                                  ) extends Metric

  /**
    * Represents a node in a graph of CBCs.  While building the graph a CBC may have 0-many children
    * and 0-many parents.  At the end of graph construction, if the node has many parents that lead
    * back to > 1 single root, the set of parents is cleared.
    *
    * @param cb the cell barcode at this node
    * @param count the count of UMIs observed by this cell barcode
    */
  final case class Node(cb: Array[Byte], count: Int) {
    val parents: mutable.Buffer[Node]  = new mutable.ArrayBuffer[Node](2)
    val children: mutable.Buffer[Node] = new mutable.ArrayBuffer[Node](2)

    /** Total count of UMIs including child counts. Note may double count where a child observed the same UMI. */
    def total: Int = count + children.map(_.total).sum

    /** The string barcode of the CBC. */
    def barcode: String = new String(cb)

    /** The set of roots for this node. The set will never be empty; if there are no parents the node is it's own root. */
    def roots: Set[Node] = {
      if (parents.isEmpty) Set(this)
      else parents.toSet.flatMap { p: Node => p.roots }
    }
  }
}

@clp(group=ClpGroups.Personal, description=
  """
    |QCs cell barcodes and UMIs in DropSeq data.
  """)
class CheckDropSeqData
(@arg(flag='i', doc="Input fastq file containing the cell barcode and UMI.") val input: PathToFastq,
 @arg(flag='o', doc="Output file prefix.") val output: PathPrefix,
 @arg(flag='c', doc="Length on the cell barcode.") val cellBarcodeLength: Int = 12,
 @arg(flag='u', doc="Length of the UMI.") val umiLength: Int,
 @arg(flag='C', doc="Expected cell count.") val cellCount: Int,
 @arg(flag='d', doc="Description of sample/experiment.") val description: String,
 @arg(flag='q', doc="Minimum acceptable base quality.") val minBaseQuality: Int = 10,
 @arg(flag='f', doc="Fold difference between parent and child barcodes to combine.") val fold: Double = 4.0,
 @arg(flag='g', doc="Graph-correct the cell barcodes and create corrected outputs.") val graphCorrect: Boolean = true,
 @arg(flag='n', doc="Generate detailed information on the top N cell barcodes.") val topN: Option[Int] = None,
 @arg(flag='t', doc="Threads to use when comparing barcodes.") val threads: Int = 8
) extends FgBioTool with LazyLogging {
  import CheckDropSeqData._


  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in = FastqSource(input)
    val iter = in.iterator.bufferBetter
    val readLen = iter.headOption.map(_.length).getOrElse(cellBarcodeLength + umiLength)

    val usableCbLength  = min(cellBarcodeLength, readLen)
    val usableUmiLength = max(0, min(umiLength, readLen-usableCbLength))
    val totalLength     = usableCbLength + usableUmiLength

    require(usableCbLength > 0 && usableUmiLength > 0,
      s"Read length is not long enough. Usable CB and UMI lengths: $usableCbLength + $usableUmiLength")

    if (usableCbLength != cellBarcodeLength || usableUmiLength != umiLength) {
      logger.warning(s"Specified cell barcode length $cellBarcodeLength and umi length $umiLength.")
      logger.warning(s"Total read length is $readLen. Using CB Length=$usableCbLength and UMI Length=$usableUmiLength.")
    }

    ///////////////////////////////////////////////////////////////////////////
    // Read the data in and make counts of unique UMIs per cell barcode
    ///////////////////////////////////////////////////////////////////////////
    var total = 0L
    val counter = new SimpleCounter[String]
    val progress = ProgressLogger(logger, unit=5e6.toInt)

    iter.foreach { rec =>
      total += 1
      val bases = rec.bases.substring(0, totalLength)
      val quals = QualityEncoding.Standard.toStandardNumeric(rec.quals.substring(0, totalLength))

      if (bases.forall(c => c != 'N' && c != '.') && quals.forall(_ >= minBaseQuality)) {
        val cb  = bases.substring(0, usableCbLength)
        val umi = bases.substring(usableCbLength, usableCbLength + usableUmiLength)
        counter.count(cb + "-" + umi)
      }

      progress.record()
    }

    in.close()
    logger.info(f"Read ${total}%,d reads. Found ${counter.total}%,d passing reads.")

    ///////////////////////////////////////////////////////////////////////////
    // Now do something interesting with the data
    ///////////////////////////////////////////////////////////////////////////
    val metrics = new ArrayBuffer[DropSeqSummaryMetrics]()
    val (umiCbcCounter, rawMetrics) = summarize(counter, "raw", s => s)
    metrics += rawMetrics
    writeOutputs(umiCbcCounter, output + ".raw", "raw")

    val nodes = umiCbcCounter.map { case (cb, n) => Node(cb.getBytes(), n.toInt) }.toIndexedSeq.sortBy(- _.count)
    val nodeSubset = if (nodes.size < 100000) nodes else nodes.filter(_.count >= 5)

    if (graphCorrect) {
      buildGraph(nodeSubset, this.fold)
      val map = nodes.iterator.map(n => n.barcode -> n.roots.head.barcode).toMap
      val (correctedCounter, correctedMetrics) = summarize(counter, "corrected", map.apply)
      metrics += correctedMetrics
      writeOutputs(correctedCounter, output + ".corrected", "corrected")
    }

    Metric.write(Paths.get(output + ".summary_metrics.txt"), metrics)

    topN.foreach { n =>
      nodes.foreach { node => node.parents.clear(); node.children.clear() }
      generateDebugInfo(nodes=nodeSubset, topCbcs=n, prefix=Paths.get(output + ".debug"))
    }
  }


  /**
    * Generates summary data from a counter that has counted CBC-UMI combinations.
    *
    * @param counter a counter of CBC-UMI -> count of reads
    * @param dataType a short description of the data type (e.g. raw, corrected)
    * @param mapping a mapping from raw cell barcode to corrected cell barcode
    * @return a tuple containing a counter of CBC->UMIs and a set of summary metrics
    */
  def summarize(counter: SimpleCounter[String], dataType: String,  mapping: String => String): (SimpleCounter[String], DropSeqSummaryMetrics) = {
    val readCbcCounter = new SimpleCounter[String]
    val umisByCbc      = new mutable.HashMap[String, mutable.Set[String]] with mutable.MultiMap[String,String]
    val umiCbcCounter  = new SimpleCounter[String]

    // Process the counts into counts of raw reads and UMIs per cell barcode
    counter.foreach { case (barcodes, n) =>
      val Array(cbc, umi) = barcodes.split('-')
      val correctedCbc = mapping(cbc)

      readCbcCounter.count(correctedCbc, n)
      umisByCbc.addBinding(correctedCbc, umi)
    }

    umisByCbc.foreach { case (cbc, umis) => umiCbcCounter.count(cbc, umis.size) }
    umisByCbc.clear()
    logger.info(f"${dataType}: Observed ${readCbcCounter.size}%,d cell barcodes with ${readCbcCounter.total}%,d total reads from ${umiCbcCounter.total}%,d UMIs.")

    val top10k  = umiCbcCounter.toSeq.sortBy(- _._2).take(10000).map { case (_, n) => n.toDouble }
    val top5k   = top10k.take(5000)
    val top1500 = top10k.take(1500)
    val top1000 = top10k.take(1000)
    val top500  = top10k.take(5000)

    def median(xs: Seq[Double]): Double = (xs(xs.length/2) + xs(xs.length/2 + 1)) / 2

    val metrics = DropSeqSummaryMetrics(
      sample                        = this.description,
      expected_cells                = this.cellCount,
      data_type                     = dataType,
      read_count                    = counter.total,
      cbc_count                     = umiCbcCounter.size,
      umi_count                     = umiCbcCounter.total,
      median_umis_top_500_cells     = median(top500),
      mean_umis_top_500_cells       = top500.sum / top500.size,
      median_umis_top_1000_cells    = median(top1000),
      mean_umis_top_1000_cells      = top1000.sum / top1000.size,
      median_umis_top_1500_cells    = median(top1500),
      mean_umis_top_1500_cells      = top1500.sum / top1500.size,
      median_umis_top_5000_cells    = median(top5k),
      mean_umis_top_5000_cells      = top5k.sum / top5k.size,
      median_umis_top_10000_cells   = median(top10k),
      mean_umis_top_10000_cells     = top10k.sum / top10k.size,
      fraction_cbcs_with_250_umis   = umiCbcCounter.count(_._2 >=  250) / umiCbcCounter.size.toDouble,
      fraction_cbcs_with_500_umis   = umiCbcCounter.count(_._2 >=  500) / umiCbcCounter.size.toDouble,
      fraction_cbcs_with_1000_umis  = umiCbcCounter.count(_._2 >= 1000) / umiCbcCounter.size.toDouble,
    )

    (umiCbcCounter, metrics)
  }

  /**
    * Writes out a detailed file giving UMI counts for every CBC, and then generates a plot.
    */
  def writeOutputs(umiCbcCounter: SimpleCounter[String], path: String, dataType: String): Unit = {
    val countsPath = Paths.get(path + ".counts.txt")
    val plotPath   = Paths.get(path + ".counts.png")
    val out = new PrintStream(countsPath.toFile)
    out.println("cell_barcode\tdata_type\tcount")

    umiCbcCounter.toSeq.sortBy(- _._2).foreach { case (cbc, umis) =>
      out.println(s"${cbc}\t${dataType}\t${umis}")
    }

    out.close()
    Rscript.execIfAvailable(PlotScript, countsPath.toString, plotPath.toString, s"${description} ${dataType}" , cellCount.toString)
  }


  /**
    * Builds a directed adjacency graph between the CBCs represented by `nodes`.  A link is made from
    * `parent` to `child` if `parent.count >= child.count * fold` AND the two barcodes are at most
    * a single mismatch different.
    *
    * After building the graph, any node that leads back to multiple roots is orphaned (i.e. has all
    * links to parents broken) as we cannot tell which original CBC originated the child.
    *
    * @param nodes the set of nodes to build the graph between
    * @param fold the minimum fold-difference between parent and child counts
    */
  def buildGraph(nodes: IndexedSeq[Node], fold: Double): Unit = {
    logger.info(f"Building graph from ${nodes.size}%,d nodes.")
    var connectionsMade = new AtomicLong(0)
    var multiParentSingleRoots = 0L
    var multiParentMultiRoots  = 0L
    val progress = ProgressLogger(logger, noun="nodes", unit=5000)

    // Build a graph where each pair of nodes is connected if they are within a small edit distance
    // and the child has `fold` less UMIs than the parent
    nodes.zipWithIndex.parWith(threads).foreach { case (parent, index) =>
      val maxChildCount = parent.count / fold

      forloop(from=index+1, until=nodes.length) { j =>
        val child = nodes(j)
        if (child.count <= maxChildCount && closeEnoughMismatchOnly(parent.cb, child.cb)) {
          parent.synchronized { parent.children += child }
          child.synchronized  { child.parents += parent }
          connectionsMade.incrementAndGet()
        }
      }

      progress.synchronized(progress.record())
    }

    // Now go through and break connections where a node has multiple parents
    nodes.iterator.filter(_.parents.size > 1).foreach { node =>
      val roots = node.roots
      if (roots.size == 1) {
        val best = node.parents.maxBy(_.count)
        node.parents.filter(_ != best).foreach(_.children -= node)
        node.parents.clear()
        node.parents += best
        multiParentSingleRoots += 1
      }
      else {
        node.parents.foreach(p => p.children -= node)
        node.parents.clear()
        multiParentMultiRoots += 1
      }
    }

    logger.info(f"Made ${connectionsMade.get()}%,d connections.")
    logger.info(f"${multiParentSingleRoots}%,d nodes has multiple parents but a single root.")
    logger.info(f"${multiParentMultiRoots}%,d nodes has multiple parents and multiple roots (and were orphaned).")
  }


  /**
    * Builds an alternative directed adjacency graph that uses a gapped alignment between cell barcodes, allowing
    * for up to two mismatches and/or inserted/deleted bases.  Because this method is much slower it is limited
    * to processing the top `topCbcs` plus their children.  However, even a DAG built from a small number of nodes
    * as a starting point can expand rapidly it is also limited to finding parents for the `maxTotal` nodes
    * during traversal.
    *
    * The result is an output file that contains trees/subgraphs for the `topCbcs` for exploratory purposes.
    *
    * @param nodes the set of nodes to build graphs between
    * @param topCbcs the number of top CBCs to start graphs from
    * @param maxTotal the maximum number of CBCs to try and find the children of
    * @param prefix the path prefix of the output file
    */
  def generateDebugInfo(nodes: IndexedSeq[Node],
                        topCbcs: Int,
                        maxTotal: Int = 10000,
                        prefix: PathPrefix): Unit = {
    logger.info(s"Building debug graph for top $topCbcs nodes and children.")
    val aligner = Aligner(4, -4, 0, -2, mode=Global)
    val progress = ProgressLogger(logger, noun="debug graph nodes", unit=20)

    // Build a graph where each pair of nodes is connected if they are within a small edit distance
    // and the child has `fold` less UMIs than the parent
    val toExamine = new java.util.concurrent.ConcurrentLinkedDeque[Node]()
    nodes.take(topCbcs).foreach(toExamine.add)
    var processed = 0

    while (toExamine.nonEmpty && processed <= maxTotal) {
      val ns = new mutable.ArrayBuffer[Node]
      Range(0, threads*3).foreach { _ => if (!toExamine.isEmpty) ns += toExamine.removeFirst() }
      processed += ns.size

      ns.filter(_.children.isEmpty).parWith(threads).foreach { parent =>
        val maxChildCount = (parent.count / 2.0).toInt
        val index = nodes.indexWhere(_ eq parent)

        forloop(from=index+1, until=nodes.length) { j =>
          val child = nodes(j)
          if (child.parents.isEmpty && child.count <= maxChildCount && closeEnoughGappedAlignment(parent.cb, child.cb, aligner)) {
            parent.synchronized { parent.children += child }
            child.synchronized  { child.parents += parent }
            toExamine.add(child)
          }
        }

        progress.record()
      }
    }

    val out = new PrintStream(prefix + ".txt")

    def print(node: Node, parent: Option[Node] = None, depth: Int = 0): Unit = {
      val prefix = " " * (4 * depth)
      val cig = parent match {
        case None    => ""
        case Some(p) => ". Cigar:" + aligner.align(node.barcode, p.barcode).cigar.toString()
      }

      val otherParents = node.parents.filterNot(p => parent.contains(p)).map(_.barcode)
      val otherString = if (otherParents.isEmpty) "" else ". Other parents: " + otherParents.mkString(",")

      out.println(s"$prefix${node.barcode} ${node.count}/${node.total}${cig}${otherString}")
      node.children.foreach(c => print(c, parent=Some(node), depth=depth+1))
    }

    nodes.take(topCbcs).foreach(print(_))
    out.close()
  }

  /** Returns true if the two sequences are at most one mismatch apart, else false. */
  @inline final def closeEnoughMismatchOnly(s1: Array[Byte], s2: Array[Byte]): Boolean = {
    var mismatches = 0
    val len = s1.length
    var i = 0

    while (i < len && mismatches < 2) {
      if (s1(i) != s2(i)) mismatches += 1
      i += 1
    }

    mismatches <= 1
  }

  /**
    * Returns true if the two sequences have at least `length-2` matching bases and:
    *   - There is at most one insertion/deleteion (and matching deletion/insertion due to global alignment)
    *   - Either the opening indel is at the start or the closing indel is at the end of the alignment
    */
  final def closeEnoughGappedAlignment(s1: Array[Byte], s2: Array[Byte], aligner: Aligner): Boolean = {
    val alignment = aligner.align(s1, s2)
    val cig = alignment.cigar
    val matching = cig.iterator.filter(_.operator == CigarOperator.EQ).sumBy(_.length)
    val indels = cig.count(c => c.operator == CigarOperator.I || c.operator == CigarOperator.D)
    val firstIsIndel = indels > 0 && (cig.head.operator == CigarOperator.I || cig.head.operator == CigarOperator.D)
    val lastIsIndel  = indels > 0 && (cig.last.operator == CigarOperator.I || cig.last.operator == CigarOperator.D)
    val indelsOk = (indels == 0) || (indels == 2 && (firstIsIndel || lastIsIndel))

    matching >= s1.length - 2 && indelsOk
  }
}
