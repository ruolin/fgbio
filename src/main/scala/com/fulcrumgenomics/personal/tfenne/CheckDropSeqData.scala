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
import com.fulcrumgenomics.util.{Io, ProgressLogger, Rscript}
import htsjdk.samtools.CigarOperator

import math.{max, min}
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer

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

  private val PlotScript = "com/fulcrumgenomics/personal/tfenne/checkDropSeqData.R"

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  case class DropSeqSummaryMetrics(sample: String,
                                   expected_cells: Long,
                                   data_type: String,
                                   read_count: Long,
                                   cbc_count: Long,
                                   umi_count: Long,
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
                                  )


  final case class Node(cb: Array[Byte], count: Int) {
    val parents: mutable.Buffer[Node]  = new mutable.ArrayBuffer[Node](2)
    val children: mutable.Buffer[Node] = new mutable.ArrayBuffer[Node](2)

    def total: Int = count + children.map(_.total).sum
    def barcode: String = new String(cb)
    def roots: Set[Node] = {
      if (parents.isEmpty) Set(this)
      else parents.toSet.flatMap { p: Node => p.roots }
    }
  }

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
    var passing = 0L
    val rawCellCounter = new SimpleCounter[String]
    val uqCellCounter  = new SimpleCounter[String]
    val sequences = new mutable.HashSet[String]
    val progress = ProgressLogger(logger, unit=5e6.toInt)


    iter.foreach { rec =>
      total += 1
      val bases = rec.bases.substring(0, totalLength)
      val quals = QualityEncoding.Standard.toStandardNumeric(rec.quals.substring(0, totalLength))

      if (bases.forall(c => c != 'N' && c != '.') && quals.forall(_ >= minBaseQuality)) {
        passing += 1
        val cb  = bases.substring(0, usableCbLength)
        val umi = bases.substring(usableCbLength, usableCbLength + usableUmiLength)

        val novel = sequences.add(cb + umi)
        rawCellCounter.count(cb)
        if (novel) uqCellCounter.count(cb)
      }

      progress.record()
    }

    in.close()

    ///////////////////////////////////////////////////////////////////////////
    // Now do something interesting with the data
    ///////////////////////////////////////////////////////////////////////////

    def printInfo(nodes: Seq[Node], prefix: String): Unit = {
      val roots = nodes.filter(_.parents.isEmpty)
      logger.info(f"$prefix cell barcodes with UMIs >=0, 5, 25: ${roots.size}%,d, ${roots.count(_.total >= 5)}%,d, ${roots.count(_.total >= 25)}%,d")
    }

    logger.info(f"Read ${total}%,d reads. Found ${passing}%,d passing reads.")
    logger.info(f"Observed ${rawCellCounter.size}%,d cell barcodes with ${rawCellCounter.total}%,d total reads from ${uqCellCounter.total}%,d UMIs.")

    val nodes = uqCellCounter.map { case (cb, n) => Node(cb.toUpperCase().getBytes(), n.toInt) }.toIndexedSeq.sortBy(- _.count)
    val nodeSubset = if (nodes.size < 100000) nodes else nodes.filter(_.count >= 5)
    printInfo(nodes, "Raw")

    if (graphCorrect) {
      writeOutputs(nodes, output + ".raw", description + " Raw", cellCount)
      buildGraph(nodeSubset, this.fold)
      printInfo(nodes, "Corrected")
      writeOutputs(nodes, output + ".corrected", description + " Corrected", cellCount)
    }

    topN.foreach { n =>
      nodes.foreach { node => node.parents.clear(); node.children.clear() }
      generateDebugInfo(nodeSubset, n, Paths.get(output + ".debug"))
    }
  }


  /**
    *
    * @param nodes
    * @param path
    * @param description
    * @param cellCount
    */
  def writeOutputs(nodes: Seq[Node], path: String, description: String, cellCount: Int): Unit = {
    val countsPath = Paths.get(path + ".counts.txt")
    val plotPath   = Paths.get(path + ".counts.png")
    val out = new PrintStream(countsPath.toFile)
    out.println("cell_barcode\tcount\texact_count")

    nodes.filter(_.parents.isEmpty).sortBy(-_.total).foreach { node =>
      out.println(s"${node.barcode}\t${node.total}\t${node.count}")
    }

    out.close()
    Rscript.execIfAvailable(PlotScript, countsPath.toString, plotPath.toString, description, cellCount.toString)
  }


  /**
    *
    * @param nodes
    * @param fold
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
    logger.info(f"${multiParentMultiRoots}%,d nodes has multiple parents and multiple roots (and were pruned)..")
  }


  /**
    *
    * @param nodes
    * @param howMany
    * @param prefix
    */
  def generateDebugInfo(nodes: IndexedSeq[Node], howMany: Int, prefix: PathPrefix): Unit = {
    logger.info(s"Building debug graph for top $howMany nodes and children.")
    val aligner = Aligner(4, -4, 0, -2, mode=Global)
    val progress = ProgressLogger(logger, noun="debug graph nodes", unit=20)

    // Build a graph where each pair of nodes is connected if they are within a small edit distance
    // and the child has `fold` less UMIs than the parent
    val toExamine = new java.util.concurrent.ConcurrentLinkedDeque[Node]()
    nodes.take(howMany).foreach(toExamine.add)
    var processed = 0

    while (toExamine.nonEmpty && processed <= 50000) {
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

    nodes.take(howMany).foreach(print(_))
    out.close()
  }

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
