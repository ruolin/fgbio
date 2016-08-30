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

import java.io.PrintStream
import java.net.URL
import java.nio.file.Paths
import java.text.SimpleDateFormat
import java.util.Date

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.metagenomics.NcbiTaxonomy.TaxonId
import com.fulcrumgenomics.util.{DelimitedDataParser, Io}
import dagr.commons.io.PathUtil
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}

import scala.collection.mutable.ListBuffer
import scala.collection.parallel.ForkJoinTaskSupport
import scala.concurrent.forkjoin.ForkJoinPool
import scala.io.Source

/** Enumeration of assembly levels in order of goodness. */
private object AssemblyLevel extends Enumeration {
  val CompleteGenome, Chromosome, Scaffold, Contig = Value
  def apply(s: String): Value = withName(s.replace(" ", ""))
}

/** Enumeration of Genome Representation levels in order of goodness. */
private object GenomeRepresentation extends Enumeration {
  val Full, Partial = Value
  def apply(s: String): Value = withName(s)
}

/** Enumeration of RefseqCategories in order of goodness. */
private object RefseqCategory extends Enumeration {
  val reference_genome, representative_genome, na = Value
  def apply(s: String): Value = withName(s.replace(' ', '_'))
}

/** Case class representing an NCBI registered assembly. */
private case class Assembly(accession: String, taxId: Int, speciesTaxId: Int, organismName: String,
                            refseqCategory: RefseqCategory.Value, genomeRepresentation: GenomeRepresentation.Value,
                            versionStatus: String, assemblyLevel: AssemblyLevel.Value,
                            releaseDate: Date, ftpPath: String)

/** Case class representing the taxons we want representation of. */
private case class DesiredTaxon(name: String, id: TaxonId, assemblies: ListBuffer[Assembly] = ListBuffer())

@clp(group = ClpGroups.Metagenomics, description =
  """
    |Finds representative genomes for taxa from the NCBI assembly database
    |and downloads them.
  """)
class FindGenomes
( @arg(flag="i", doc="Tab-separated file with two columns, 'name' and 'taxid'.") val input: FilePath,
  @arg(flag="o", doc="Output directory to write to.") val output: DirPath = Paths.get(System.getProperty("user.dir")),
  @arg(flag="r", doc="Use RefSeq (default is more inclusive GenBank).") val refseq: Boolean = false,
  @arg(flag="t", doc="Directory with NCBI taxonomy files (names.dmp and nodes.dmp).") val taxonomyDir: DirPath,
  @arg(flag="d", doc="Download chosen genomes after writing reports.") val download: Boolean = false,
  @arg(flag="m", doc="When downloading genomes, hard-mask soft-masked sequence.") val hardMask: Boolean = true,
  @arg(flag="T", doc="Threads to use when downloading genomes.") val threads: Int = 1
) extends FgBioTool with LazyLogging {

  // Check the inputs
  Io.assertReadable(input)
  Io.assertWritableDirectory(output)
  private val Seq(names, nodes) = Seq("names.dmp", "nodes.dmp").map(taxonomyDir.resolve)
  Io.assertReadable(Seq(names, nodes))

  private val BaseUrl = "ftp://ftp.ncbi.nlm.nih.gov/genomes/" + (if (refseq) "refseq" else "genbank")
  private val SummaryFilename = "assembly_summary.txt"

  /* Sub-directories of the NCBI genbank and refseq assembly FTPs */
  private val SubDirs = Seq("archaea", "bacteria", "fungi", "invertebrate", "metagenomes", "other",
                            "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other")

  override def execute(): Unit = {
    val taxonsById = loadInput
    val assemblies = loadSummaries
    logger.info("Loading NCBI taxonomy files.")
    val taxdb = NcbiTaxonomy(names, nodes)

    logger.info("Locating assemblies for desired taxons.")
    populateAssemblies(taxonsById, taxdb, assemblies)

    taxonsById.values.foreach { dt =>
      logger.info(s"Found ${dt.assemblies.size} assemblies for ${dt.name} (${dt.id}).")
    }

    val sortedTaxons = taxonsById.values.toSeq.sortBy(_.name)
    outputSummaries(sortedTaxons, taxdb)
    if (download) fetchGenomes(taxons=sortedTaxons, dir=output, hardMask=hardMask)
  }

  /** Loads up the set of desired taxons from the input file. */
  private def loadInput: Map[TaxonId,DesiredTaxon] = {
    DelimitedDataParser(input, '\t')
      .map(row => DesiredTaxon(name=row[String]("name"), id=row[Int]("taxid")))
      .map(t => t.id -> t)
      .toMap
  }

  /**
    * Fetches all the summary files from NCBI and loads up the list of assemblies for us
    * to query.
    * @return a Seq of assembly objects
    */
  private def loadSummaries: Seq[Assembly] = {
    logger.info("Loading assembly summary files from NCBI.")
    val dateFmt = new SimpleDateFormat("yyyy/MM/dd")
    val results = ListBuffer[Assembly]()

    SubDirs.foreach { dir =>
      val url = new URL(BaseUrl + "/" + dir + "/" + SummaryFilename)
      // The files have a comment line, then a header line which starts with "# ", so we need to strip that out
      val lines = Source.fromURL(url).getLines().toList match {
        case junk :: headers :: rest => headers.substring(2) :: rest
        case xs =>
          logger.warning(s"${SummaryFilename} for ${dir} contained less than two rows.")
          Nil
      }

      val assemblies = DelimitedDataParser(lines, '\t').map(row => Assembly(
        accession            = row[String]("assembly_accession"),
        taxId                = row[Int]("taxid"),
        speciesTaxId         = row[Int]("species_taxid"),
        organismName         = row[String]("organism_name"),
        refseqCategory       = RefseqCategory(row[String]("refseq_category")),
        genomeRepresentation = GenomeRepresentation(row[String]("genome_rep")),
        versionStatus        = row[String]("version_status"),
        assemblyLevel        = AssemblyLevel(row[String]("assembly_level").replace(" ", "")),
        releaseDate          = dateFmt.parse(row.get[String]("seq_rel_date").getOrElse("1990/01/01")),
        ftpPath              = row[String]("ftp_path")
      )).filter(_.versionStatus == "latest").toSeq

      logger.info(s"Loaded ${assemblies.length} assemblies from ${dir}")
      results.appendAll(assemblies)
    }

    results
  }

  /**
    * Populates the assemblies value on each DesiredTaxon with all assemblies of that taxon or
    * taxons which are descendents of that taxon.
    */
  def populateAssemblies(targets: Map[TaxonId, DesiredTaxon], taxonomy: NcbiTaxonomy, assemblies: Seq[Assembly]): Unit = {
    assemblies.foreach { assembly =>
      taxonomy.get(assembly.taxId) match {
        case None => logger.warning(s"Taxid ${assembly.taxId} for organism ${assembly.organismName} and assembly ${assembly.accession} not found in Taxonomy files.")
        case Some(taxon) =>
          (Iterator(taxon) ++ taxon.ancestors) foreach { t =>
            targets.get(t.id).foreach(dt => dt.assemblies.append(assembly))
          }
      }
    }

    targets.values.foreach { dt =>
      val sorted = dt.assemblies.sortBy(a => (a.genomeRepresentation, a.refseqCategory, a.assemblyLevel, a.taxId, -a.releaseDate.getTime))
      dt.assemblies.clear()
      dt.assemblies.appendAll(sorted)
    }
  }

  /** Outputs per-taxon and per-assembly summary files. */
  def outputSummaries(targets: Seq[DesiredTaxon], taxonomy: NcbiTaxonomy): Unit = {
    def line(fields: Any*): String = fields.mkString("\t");
    def ancestors(id: TaxonId): String = taxonomy(id).ancestors.toSeq.reverse.map(_.name).mkString("->")

    {
      val out = new PrintStream(Io.toOutputStream(output.resolve("taxons.txt")))
      out.println(line("name", "requested_taxid", "requested_rank", "assemblies", "assembly", "organism", "taxid", "ancestors"))

      targets.foreach { t =>
        val a = t.assemblies.headOption
        val reqRank = taxonomy(t.id).rank.getOrElse("no rank")
        val accession = a.map(_.accession).getOrElse("no assembly")
        val organism = a.map(_.organismName).getOrElse("no assembly")
        val taxid = a.map(_.taxId).getOrElse("no assembly")
        val ancestor = a.map(x => ancestors(x.taxId)).getOrElse("no assembly")
        out.println(line(t.name, t.id, reqRank, t.assemblies.size, accession, organism, taxid, ancestor))
      }

      out.close()
    }

    {
      val fmt = new SimpleDateFormat("yyyy/MM/dd")
      val out = new PrintStream(Io.toOutputStream(output.resolve("assemblies.txt")))
      out.println(line("name", "requested_taxid", "requested_rank", "index", "assembly", "organism", "taxid", "assembly_level",
        "refseq_category", "genome_representation", "release_date", "ftp_path"))

      for (t <- targets; (a, idx) <- t.assemblies.zipWithIndex) {
        val reqRank = taxonomy(t.id).rank.getOrElse("no rank")
        val taxon = taxonomy(a.taxId)
        out.println(line(t.name, t.id, reqRank, idx + 1, a.accession, a.organismName, a.taxId, a.assemblyLevel, a.refseqCategory,
          a.genomeRepresentation, fmt.format(a.releaseDate), a.ftpPath, ancestors(a.taxId)))
      }

      out.close()
    }
  }

  /** Fetches genomes for all taxa that have genomes allocated. */
  def fetchGenomes(taxons: Iterable[DesiredTaxon], dir: DirPath, hardMask: Boolean = true): Unit = {
    val par = taxons.toSeq.par
    par.tasksupport = new ForkJoinTaskSupport(new ForkJoinPool(this.threads))

    par.foreach { taxon =>
      taxon.assemblies.toList match {
        case first :: rest =>
          val name = PathUtil.basename(first.ftpPath, trimExt=false)
          val url  = new URL(first.ftpPath + "/" + name + "_genomic.fna.gz")
          val out  = dir.resolve(first.organismName.replace(' ', '_') + "." + first.accession + ".fa")
          fetchGenome(url, out, first.taxId, hardMask=hardMask)
        case Nil => Unit
      }
    }
  }

  /** Fetches a genome from the provided URL and writes it to the target file. */
  def fetchGenome(url: URL, fasta: FilePath, taxId: TaxonId, hardMask: Boolean = true): Unit = {
    logger.info(s"Downloading assembly from: ${url}")
    val out = Io.toWriter(fasta)
    Io.toSource(url).getLines()
      .map { line =>
        if (line.startsWith(">")) {
          val fields = line.trim.split(' ')
          fields(0) + "|kraken:taxid|" + taxId + fields.iterator.drop(1).mkString(" ", " ", "")
        }
        else if (hardMask) line.map(ch => if (Character.isLowerCase(ch)) 'N' else ch)
        else line
      }
      .foreach(line => out.append(line).append('\n'))

    out.close()
  }
}
