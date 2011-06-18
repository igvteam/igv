/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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

    private Strand strand;
    private int startPosition;
    private List<AminoAcid> sequence;
    boolean nonNullSequence;

    public AminoAcidSequence(Strand strand, int startPosition, List<AminoAcid> sequence) {
        this.strand = strand;
        this.startPosition = startPosition;
        this.sequence = sequence;

        // Look for a non null sequence.  Sequences are null if the sequence
        // directory is undefined or unreachable.  
        nonNullSequence = false;
        for (AminoAcid aa : sequence) {
            if (aa != null) {
                nonNullSequence = true;
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
}
