/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.feature;

import org.broad.igv.feature.genome.Genome;

/**
 * @author Jim Robinson
 * @date 12/26/11
 */
class Codon {

    // Position in protein coordinates, first amino acid is numbered 1
    private int proteinPosition;

    // Triplet -- genomic (base) positions
    // in "zero" based (UCSC style) coordinates
    private int[] genomePositions = new int[]{-1, -1, -1};
    private int nextPos = 0;
    private AminoAcid aminoAcid;
    private Strand strand;
    private int incr = 1;
    private String chr;
    private String sequence;

    Codon(String chr, int proteinPosition) {
        this(chr, proteinPosition, Strand.POSITIVE);
    }

    public Codon(String chr, int proteinPosition, Strand strand) {
        this.proteinPosition = proteinPosition;
        this.strand = strand;
        this.chr = chr;
        if (this.strand == Strand.NEGATIVE) {
            this.nextPos = 2;
            this.incr = -1;
        }
    }

    /**
     * If on positive strand, will add each successive
     * gp at end. If on negative strand, they will be added in reverse order.
     *
     * @param gp
     */
    void setNextGenomePosition(int gp) {
        try {
            genomePositions[nextPos] = gp;
            nextPos += incr;
        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("Genome positions have already been set, cannot set any more");
        }
    }

    boolean isGenomePositionsSet() {
        boolean set = true;
        for (int ii : genomePositions) {
            set &= ii >= 0;
        }
        return set;
    }

    public void setAminoAcid(AminoAcid aa) {
        this.aminoAcid = aa;
    }

    public int getProteinPosition() {
        return proteinPosition;
    }

    public int[] getGenomePositions() {
        return genomePositions;
    }

    public AminoAcid getAminoAcid() {
        return aminoAcid;
    }

    public void calcSequence(Genome genome) {
        if (!this.isGenomePositionsSet()) {
            throw new IllegalStateException("Must set genome positions first");
        }
        int[] positions = this.getGenomePositions();
        String aas = "";
        for (int start : positions) {
            final byte[] nucSequence = genome.getSequence(chr, start, start + 1);
            if (nucSequence == null) {
                // No sequence.
            } else {
                aas += new String(nucSequence);
            }
        }

        if (strand == Strand.NEGATIVE) {
            aas = AminoAcidManager.getNucleotideComplement(aas);
        }

        this.sequence = aas;
    }

    String getSequence() {
        return sequence;
    }
}
