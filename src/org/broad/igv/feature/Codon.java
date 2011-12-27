package org.broad.igv.feature;

/**
 * @author Jim Robinson
 * @date 12/26/11
 */
class Codon {

    private int proteinPosition;   // Position in protein coordinates, first amino acid is numbered 1
    private int genomePosition1 = -1;  // Triplet -- genomic (base) positions
    private int genomePosition2 = -1;  //    in "zero" based (UCSC style) coordinates
    private int genomePosition3 = -1;
    private AminoAcid aminoAcid;

    Codon(int proteinPosition) {
        this.proteinPosition = proteinPosition;
    }

    void setNextGenomePosition(int gp) {
        if (genomePosition1 < 0) {
            genomePosition1 = gp;
        } else if (genomePosition2 < -0) {
            genomePosition2 = gp;
        } else if (genomePosition3 < 0) {
            genomePosition3 = gp;
        }
        else {
            // TODO -- throw an exception?
        }
    }

    boolean isGenomePositionsSet() {
        return genomePosition1 >= 0 && genomePosition2 >= 0 && genomePosition3 >= 0;
    }

    public void setAminoAcid(AminoAcid aa) {
        this.aminoAcid = aa;
    }

    public int getProteinPosition() {
        return proteinPosition;
    }

    public int getGenomePosition1() {
        return genomePosition1;
    }

    public int getGenomePosition2() {
        return genomePosition2;
    }

    public int getGenomePosition3() {
        return genomePosition3;
    }

    public AminoAcid getAminoAcid() {
        return aminoAcid;
    }
}
