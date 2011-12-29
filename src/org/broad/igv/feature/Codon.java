package org.broad.igv.feature;

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

    Codon(int proteinPosition) {
        this.proteinPosition = proteinPosition;
    }

    void setNextGenomePosition(int gp) {
        if (nextPos < 3) {
            genomePositions[nextPos] = gp;
            nextPos++;
        } else {
            throw new IllegalStateException("Genome positions have already been set, cannot set any more");
        }
    }

    boolean isGenomePositionsSet() {
        boolean set = true;
        for (int ii : genomePositions) {
            set &= ii > 0;
        }
        return set;
    }

    public void setAminoAcid(AminoAcid aa) {
        this.aminoAcid = aa;
    }

    public int getProteinPosition() {
        return proteinPosition;
    }

    public int getGenomePosition(int pos) {
        try {
            return genomePositions[pos];
        } catch (IndexOutOfBoundsException e) {
            throw new IllegalArgumentException("pos argument must be 0-2");
        }
    }

    public AminoAcid getAminoAcid() {
        return aminoAcid;
    }
}
