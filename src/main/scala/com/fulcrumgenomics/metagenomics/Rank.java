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

package com.fulcrumgenomics.metagenomics;

/**
 * Enum representing the ranks of taxonomy at NCBI. The value at NCBI use spaces which are replaced
 * with underscores here.  In addition the value "class" has a single underscore appended so as not
 * to conflict with the Java keyword.
 *
 * Entries are ordered such that Ranks with a lower ordinal value are always higher in the taxanomic
 * tree that values with higher orderinals!
 */
public enum Rank {
    superkingdom,
    kingdom,
    subkingdom,

    superphylum,
    phylum,
    subphylum,

    superclass,
    class_,
    subclass,
    infraclass,

    superorder,
    order,
    suborder,
    infraorder,
    parvorder,

    superfamily,
    family,
    subfamily,
    tribe,
    subtribe,

    genus,
    subgenus,

    species_group,
    species_subgroup,
    species,
    subspecies,
    varietas,
    forma
    ;

    /** Should be used in preference to name() to generate the string used by NCBI. */
    @Override public String toString() { return name().replace('_', ' ').trim(); }

    /**
     * Returns true if this rank is higher than the other rank. For example:
     *     genus.isAbove(species) == true
     *     genus.isAbove(family)  == false
     *     genus.isAbove(genus)   == false
     *
     * @see #isAtOrAbove(Rank)
     */
    public boolean isAbove(final Rank other) { return this.ordinal() < other.ordinal(); }

    /**
     * Returns true if this rank is the same as, or higher than, the other rank.
     *
     * @see #isAbove(Rank)
     */
    public boolean isAtOrAbove(final Rank other) { return this == other || isAbove(other); }

    /** Mangles the name as needed before calling valueOf(). */
    public static Rank of(final String name) {
        String mangled = name.replace(' ', '_');
        if (mangled.equals("class")) mangled += "_";
        return valueOf(mangled);
    }
}
