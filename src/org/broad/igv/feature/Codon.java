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
    private Strand strand;
    private int incr = 1;

    Codon(int proteinPosition) {
        this(proteinPosition, Strand.POSITIVE);
    }

    public Codon(int proteinPosition, Strand strand) {
        this.proteinPosition = proteinPosition;
        this.strand = strand;
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

    public int[] getGenomePositions() {
        return genomePositions;
    }

    public AminoAcid getAminoAcid() {
        return aminoAcid;
    }
}
