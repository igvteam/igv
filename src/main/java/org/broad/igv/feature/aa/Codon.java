package org.broad.igv.feature.aa;

import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.SequenceTrack;

/**
 * @author Jim Robinson
 * @date 12/26/11
 */
public class Codon {

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
    public void setNextGenomePosition(int gp) {
        try {
            genomePositions[nextPos] = gp;
            nextPos += incr;
        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("Genome positions have already been set, cannot set any more");
        }
    }

    public boolean isGenomePositionsSet() {
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
            aas = SequenceTrack.getReverseComplement(aas);
        }

        this.sequence = aas;
    }

    public String getSequence() {
        return sequence;
    }
}
