/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author jrobinso
 */
public class AminoAcidManager {

    static String filename = "codonTable.txt";

    static Map<String, AminoAcid> codonTable;

    /**
     * Get the amino acid sequence for an interval.
     * Assumptions and conventions
     * <p/>
     * The start and end positions are on the positive strand (reference coordinate system)
     * irrespective of the read direction.
     * <p/>
     * Reading will begin from the startPosition if strand == POSITIVE, endPosition if NEGATIVE
     *
     * @param startPosition
     * @param endPosition
     * @return
     */
    public static AminoAcidSequence getAminoAcidSequence(String genomeId, String chr,
                                                         int startPosition, int endPosition, Strand strand) {

        byte[] seqBytes = SequenceManager.readSequence(genomeId, chr, startPosition, endPosition);
        if (seqBytes == null) {
            return null;
        } else {
            String nucSequence = new String(seqBytes);
            List<AminoAcid> acids = getAminoAcids(nucSequence, strand);

            // Set the start position of this amino acid.  
            int aminoStart = startPosition;
            return new AminoAcidSequence(strand, aminoStart, acids);
        }

    }

    /**
     * Return an amino acid sequence for the input sequence.
     *
     * @param sequence
     * @param direction
     * @return
     */
    public static List<AminoAcid> getAminoAcids(String sequence, Strand direction) {

        if (codonTable == null) {
            initTable();
        }

        // Sequence must be divisible by 3. It is the responsibilty of the 
        // calling program to send a sequence properly aligned.  
        int readLength = sequence.length() / 3;
        List<AminoAcid> acids = new ArrayList(readLength);

        for (int i = 0; i <= sequence.length() - 3; i += 3) {
            String codon = sequence.substring(i, i + 3);
            if (direction == Strand.NEGATIVE) {
                codon = getNucleotideComplement(codon);
            }
            if (codonTable == null || codon == null) {
                acids.add(AminoAcid.NULL_AMINO_ACID);
            } else {
                AminoAcid aa = codonTable.get(codon);
                if (aa == null) {
                    aa = AminoAcid.NULL_AMINO_ACID;
                }
                acids.add(codonTable.get(codon));
            }

        }
        return acids;
    }

    private static String getNucleotideComplement(String sequence) {
        char[] complement = new char[sequence.length()];
        for (int i = 0; i < sequence.length(); i++) {
            int j = sequence.length() - i - 1;
            char c = sequence.charAt(i);
            switch (c) {
                case 'T':
                    complement[j] = 'A';
                    break;
                case 'A':
                    complement[j] = 'T';
                    break;
                case 'C':
                    complement[j] = 'G';
                    break;
                case 'G':
                    complement[j] = 'C';
                    break;
                default:
                    complement[j] = c;
            }
        }
        return new String(complement);
    }

    static synchronized void initTable() {

        if (codonTable == null) {


            try {

                InputStream is = AminoAcidManager.class.getResourceAsStream("/resources/" + filename);
                if (is == null) {
                    return;
                }
                codonTable = new HashMap();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                String nextLine = null;
                while ((nextLine = reader.readLine()) != null) {
                    String[] tokens = nextLine.split("\t");
                    if (tokens.length == 4) {
                        String codon = tokens[0].trim().toUpperCase();
                        String fullName = tokens[1].trim();
                        String shortName = tokens[2].trim();
                        char symbol = tokens[3].trim().charAt(0);
                        codonTable.put(codon, new AminoAcid(fullName, shortName, symbol));
                    }

                }
            } catch (IOException ex) {
                Logger.getLogger(AminoAcidManager.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
}
