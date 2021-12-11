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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.aa;

import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import com.google.gson.*;
import org.apache.logging.log4j.*;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.SequenceTrack;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * @author jrobinso
 */
public class AminoAcidManager {

    private static final Logger log = LogManager.getLogger(AminoAcidManager.class);

    /**
     * File which contains listing of amino acid names.
     * Format: Full Name \t 3 letter abbreviation \t Single letter abbrev.
     */
    private static final String AANameFilePath = "resources/AANamesTable.txt";

    /**
     * Table containing mapping from string forms (full, TLA, single-letter-abbrev)
     * to amino acid object. No codon information stored here
     */
    static final Map<String, AminoAcid> AANameMap = new HashMap<String, AminoAcid>(20);

    private static AminoAcidManager instance;

    private AminoAcidManager() {
        initAANameMap();
        try {
        } catch (JsonParseException e) {
            log.error(e);
        }
    }

    public static AminoAcidManager getInstance() {
        if (instance == null) {
                AminoAcidManager newInstance = new AminoAcidManager();
                instance = newInstance;
        }
        return instance;
    }


    /**
     * Return a list of amino acids for the input sequence of nucleotides
     *
     * @param direction
     * @param sequence
     * @return
     */

    List<CodonAA> getAminoAcids(Strand direction, String sequence, CodonTable codonTable) {

        // Sequence must be divisible by 3. It is the responsibility of the
        // calling program to send a sequence properly aligned.
        int readLength = sequence.length() / 3;
        List<CodonAA> acids = new ArrayList<>(readLength);

        if (direction == Strand.NEGATIVE) {
            sequence = SequenceTrack.getReverseComplement(sequence);
        }

        for (int i = 0; i <= sequence.length() - 3; i += 3) {
            String codon = sequence.substring(i, i + 3).toUpperCase();
            AminoAcid aa = codonTable.getAminoAcid(codon);
            CodonAA cAA = new CodonAA(codon, aa);
            acids.add(cAA);
        }

        if (direction == Strand.NEGATIVE) {
            Collections.reverse(acids);
        }

        return acids;
    }

    /**
     * Get the amino acid sequence for an interval.
     * Assumptions and conventions
     * <p>
     * The start and end positions are on the positive strand
     * irrespective of the read direction.
     * <p>
     * Reading will begin from the startPosition if strand == POSITIVE, endPosition if NEGATIVE
     *
     * @return AminoAcidSequence, or null if seqBytes == null
     */
    public AminoAcidSequence getAminoAcidSequence(Strand strand, int start, String nucSequence, CodonTable codonTable) {
        if (nucSequence == null) {
            return null;
        } else {

            int l = nucSequence.length();
            int rem = l % 3;
            int aaStart = strand == Strand.POSITIVE ? 0 : 0 + rem;

            List<CodonAA> acids = getAminoAcids(strand, nucSequence, codonTable);

            return new AminoAcidSequence(strand, start + aaStart, acids, codonTable.getId());
        }
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
        initAANameMap();

        AminoAcid aa = AANameMap.get(name);
        if (aa == null) {
            aa = AminoAcid.NULL_AMINO_ACID;
        }

        return aa;
    }

    public Set<String> getMappingSNPs(String codon, AminoAcid mutAA, CodonTable codonTable) {
        Set<String> mapSNPs = new HashSet<String>();
        Set<String> SNPs = getAllSNPs(codon);
        for (String modCodon : SNPs) {
            //We use short name because all 3 stop codon have different long names,
            //and we don't care about the difference here.
            if (codonTable.getAminoAcid(modCodon).equalsByName(mutAA.getShortName())) {
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

    /**
     * Initialize table of amino acid names, for easy lookup of
     * AminoAcid by symbols. This method is idempotent, only called once
     * to read name file.
     */
    public synchronized static void initAANameMap() {
        if (!AANameMap.isEmpty()) {
            return;
        }
        try {
            InputStream is = AminoAcidManager.class.getResourceAsStream(AANameFilePath);
            if (is == null) {
                return;
            }
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));

            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("#")) continue;
                String[] tokens = nextLine.split("\t");
                if (tokens.length == 3) {
                    String fullName = tokens[0].trim();
                    String shortName = tokens[1].trim();
                    String symbol = tokens[2].trim();
                    assert symbol.length() == 1;
                    AminoAcid aa = new AminoAcid(fullName, shortName, symbol.charAt(0));
                    for (String sym : new String[]{fullName, shortName, symbol}) {
                        if (!AANameMap.containsKey(sym)) {
                            AANameMap.put(sym, aa);
                        }
                    }
                }

            }
        } catch (IOException ex) {
            log.error(ex);
            throw new RuntimeException(ex);
        }
    }

}
