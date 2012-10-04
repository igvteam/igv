/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import java.util.List;

/**
 * Represents an amino acid sequence for an exon
 *
 * @author jrobinso
 */
public class AminoAcidSequence {

    private final Strand strand;
    private final int startPosition;
    private final List<AminoAcid> sequence;

    private boolean nonNullSequence;

    /**
     * If the sequence was computed from a nucleotide sequence,
     * we store how the calculation was performed
     */
    private final AminoAcidManager.CodonTableKey codonTableKey;

    public AminoAcidSequence(Strand strand, int startPosition, List<AminoAcid> sequence, AminoAcidManager.CodonTableKey codonTableKey) {
        this.strand = strand;
        this.startPosition = startPosition;
        this.sequence = sequence;
        this.codonTableKey = codonTableKey;

        // Look for a non null sequence.  Sequences are null if the sequence
        // directory is undefined or unreachable.  
        nonNullSequence = false;
        for (AminoAcid aa : sequence) {
            if (aa != null) {
                nonNullSequence = true;
                break;
            }
        }

    }

    public Strand getStrand() {
        return strand;
    }

    public int getStartPosition() {
        return startPosition;
    }

    public List<AminoAcid> getSequence() {
        return sequence;
    }

    public boolean hasNonNullSequence() {
        return nonNullSequence;
    }

    public AminoAcidManager.CodonTableKey getCodonTableKey() {
        return codonTableKey;
    }
}
