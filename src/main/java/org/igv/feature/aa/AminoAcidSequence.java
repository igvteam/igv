package org.igv.feature.aa;

import org.igv.feature.Strand;

import java.util.List;

/**
 * Represents an amino acid sequence for an exon
 *
 * @author jrobinso
 */
public class AminoAcidSequence {

    private final Strand strand;
    private final int start;                // Genomic position for start of sequence.
    private final List<CodonAA> sequence;

    private boolean nonNullSequence;

    /**
     * If the sequence was computed from a nucleotide sequence,
     * we store how the calculation was performed
     */
    private final Integer id;

    public AminoAcidSequence(Strand strand, int startPosition,
                             List<CodonAA> sequence,
                             Integer codonTableKey) {
        this.strand = strand;
        this.start = startPosition;
        this.sequence = sequence;
        this.id = codonTableKey;

        // Look for a non null sequence.  Sequences are null if the sequence
        // directory is undefined or unreachable.  
        nonNullSequence = false;
        for (CodonAA aa : sequence) {
            if (aa != null) {
                nonNullSequence = true;
                break;
            }
        }

    }

    public Strand getStrand() {
        return strand;
    }

    public int getStart() {
        return start;
    }

    public List<CodonAA> getSequence() {
        return sequence;
    }

    public boolean hasNonNullSequence() {
        return nonNullSequence;
    }

    public Integer getId() {
        return id;
    }

    /**
     * Return the codon sequence as a string -- primarily for debugging.
     * @return
     */
    public String getSequenceString() {
        String ss = "";
        for(CodonAA c : sequence) {
            ss += c.getAminoAcid().getSymbol();
        }
        return ss;
    }
}
