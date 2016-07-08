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

import java.nio.file.Paths

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import NcbiTaxonomyTest._

object NcbiTaxonomyTest {
  val TaxonomyDir = Paths.get("src/test/resources/com/fulcrumgenomics/metagenomics/kraken_test_db/taxonomy")
  val IdToName = Map(
    9606 -> "Homo sapiens",
    9598 -> "Pan troglodytes",
    9615 -> "Canis lupus familiaris",
    562  -> "Escherichia coli",
    5007 -> "Brettanomyces bruxellensis"
  )
}

/** Tests for the NcbiTaxonomy parsing data classes. */
class NcbiTaxonomyTest extends UnitSpec {
  def taxonomy = NcbiTaxonomy(TaxonomyDir.resolve("names.dmp"), TaxonomyDir.resolve("nodes.dmp"))

  "NcbiTaxonomy" should "load up a set of .dmp files and provide basic queries for existing values" in {
    val taxa = taxonomy
    IdToName.foreach { case (id, name) =>
      taxa(id).name shouldBe name
    }
  }

  it should "return None for IDs that do not exist in the taxonomy" in {
    val taxa = taxonomy
    (-100 to -1).foreach { id => taxa.get(id) shouldBe None }
    (Int.MaxValue-100 to Int.MaxValue).foreach { id => taxa.get(id) shouldBe None }
  }

  it should "enable finding ancestors by rank" in {
    val taxa = taxonomy
    val human = taxa(9606)
    val homo  = human.findAncestor(Rank.genus)
    homo shouldBe 'defined
    homo.get.rank.get shouldBe Rank.genus
    homo.get.name shouldBe "Homo"

    val mammalia = human.findAncestor(Rank.class_)
    mammalia shouldBe 'defined
    mammalia.get.rank.get shouldBe Rank.class_
    mammalia.get.name shouldBe "Mammalia"
  }

  it should "supporting searching for ancestors and children" in {
    val taxa = taxonomy
    val chimp = taxa(9598)
    val mammalia = chimp.ancestors.find(_.name == "Mammalia").toSeq
    mammalia.size shouldBe 1
    mammalia.head.name shouldBe "Mammalia"
    mammalia.head.children.find(_.id == 9598).size shouldBe 0
    val chimpAgain = mammalia.head.descendants.find(_.id == 9598).toSeq
    chimpAgain.length shouldBe 1
    chimpAgain.head.id shouldBe 9598
  }

  it should "correctly report whether one taxon is descended from another or not" in {
    val taxa   = taxonomy
    val chimp  = taxa(9598)
    val human  = taxa(9606)
    val mammal = taxa(40674)

    chimp.descendsFrom(mammal) shouldBe true
    chimp.descendsFrom(human) shouldBe false
    chimp.descendsFrom(chimp) shouldBe false

    human.descendsFrom(mammal) shouldBe true
    human.descendsFrom(chimp) shouldBe false
    human.descendsFrom(human) shouldBe false

    mammal.descendsFrom(mammal) shouldBe false
    mammal.descendsFrom(chimp) shouldBe false
    mammal.descendsFrom(human) shouldBe false
  }

  it should "correctly calculate level" in {
    val taxa = taxonomy
    taxa.root.level shouldBe 0
    taxa.root.children.forall(_.level == 1) shouldBe true
    taxa.root.children.forall(l1 => l1.children.forall(_.level == 2)) shouldBe true
    taxa(9606).level shouldBe taxa(9606).ancestors.length
  }

  it should "connect every node in the taxonomy tree" in {
    val taxa = taxonomy
    taxa.root.subtree.toSet           should contain theSameElementsAs taxa.taxa.toSet
    taxa.root.subtree.map(_.id).toSet should contain theSameElementsAs taxa.ids.toSet
  }

  it should "terminate ok when toString is called on a low-level taxon" in {
    taxonomy(9606).toString.contains("Homo sapiens") shouldBe true
  }
}
