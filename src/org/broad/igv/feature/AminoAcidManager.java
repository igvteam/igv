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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
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
     * The start and end positions are on the positive strand
     * irrespective of the read direction.
     * <p/>
     * Reading will begin from the startPosition if strand == POSITIVE, endPosition if NEGATIVE
     *
     * @param seqBytes
     * @param startPosition
     * @param strand
     * @return
     */
    public static AminoAcidSequence getAminoAcidSequence(byte[] seqBytes, int startPosition, Strand strand) {
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
            String codon = sequence.substring(i, i + 3).toUpperCase();
            if (direction == Strand.NEGATIVE) {
                codon = getNucleotideComplement(codon);
            }
            AminoAcid aa = getAminoAcid(codon);
            acids.add(aa);

        }
        return acids;
    }

    public static AminoAcid getAminoAcid(String codon) {
        AminoAcid aa = AminoAcid.NULL_AMINO_ACID;

        if (codonTable == null) {
            initTable();
        }

        if (codonTable != null && codon != null) {
            aa = codonTable.get(codon);
            if (aa == null) {
                aa = AminoAcid.NULL_AMINO_ACID;
            }
        }
        return aa;
    }

    /**
     * Given the 'name' of an amino acid, find a match. Lookups
     * can be by full name, short form, or single letter. Note that
     * in the case of multiple matches, the first is returned.
     * This matters most for the stop codon, whose full name
     * is ambiguous (ochre, amber, opal) if the the short form
     * or single letter is used.
     *
     * @param name
     * @return
     */
    public static AminoAcid getAminoAcidByName(String name) {
        if (codonTable == null) {
            initTable();
        }

        boolean found;
        for (AminoAcid aa : codonTable.values()) {
            found = aa.equalsByName(name);
            if (found) {
                return aa;
            }

        }

        return AminoAcid.NULL_AMINO_ACID;
    }

    public static String getNucleotideComplement(String sequence) {
        char[] complement = new char[sequence.length()];
        int jj = complement.length;
        for (int ii = 0; ii < sequence.length(); ii++) {
            char c = sequence.charAt(ii);
            jj--;
            switch (c) {
                case 'T':
                case 't':
                    complement[jj] = 'A';
                    break;
                case 'A':
                case 'a':
                    complement[jj] = 'T';
                    break;
                case 'C':
                case 'c':
                    complement[jj] = 'G';
                    break;
                case 'G':
                case 'g':
                    complement[jj] = 'C';
                    break;
                default:
                    complement[jj] = c;
            }
        }
        return new String(complement);
    }

    public static Set<String> getMappingSNPs(String codon, AminoAcid mutAA) {
        Set<String> mapSNPs = new HashSet<String>();
        Set<String> SNPs = getAllSNPs(codon);
        for (String modCodon : SNPs) {
            //todo override equals?
            //We use short name because all 3 stop codon have different long names,
            //and we don't care about the difference here.
            if (codonTable.get(modCodon).equalsByName(mutAA.getShortName())) {
                mapSNPs.add(modCodon);
            }
        }
        return mapSNPs;
    }

    /**
     * Gets all possible strings which are a SNP from
     * the provided sequence. Does not include original in
     * returned set. Assumes sequence is DNA sequence, consisting
     * of A,T,G,C, and uses that set to create SNPs.
     *
     * @param sequence
     * @return
     */
    public static Set<String> getAllSNPs(String sequence) {
        Set<String> SNPs = new HashSet<String>();
        char[] bps = "ATGC".toCharArray();
        char[] orig = sequence.toCharArray();
        char[] mod;
        for (int loc = 0; loc < orig.length; loc++) {
            mod = orig.clone();
            for (char bp : bps) {
                if (bp == orig[loc]) {
                    continue;
                }
                mod[loc] = bp;
                SNPs.add(new String(mod));
            }
        }
        return SNPs;
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
