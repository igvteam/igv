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

import org.apache.log4j.Logger;
import org.broad.igv.util.ParsingUtils;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * @author jrobinso
 */
public class AminoAcidManager {

    private static final Logger log = Logger.getLogger(AminoAcidManager.class);

    /**
     * File which contains listing of amino acid names.
     * Format: Full Name \t 3 letter abbreviation \t Single letter abbrev.
     */
    private static final String AANameFilePath = "resources/AANamesTable.txt";

    /**
     * Table containing mapping from string forms (full, TLA, single-letter-abbrev)
     * to amino acid object. No codon information stored here
     */
    private static final Map<String, AminoAcid> AANameMap = new HashMap<String, AminoAcid>(20);


    private static final String DEFAULT_CODON_TABLE_PATH = "resources/geneticCode.json";
    /**
     * File which describes how codons (e.g. ATG) map to amino acids.
     * We allow for the user loading their own, so this is not final.
     */
    private static String codonTablesPath = DEFAULT_CODON_TABLE_PATH;

    private LinkedHashMap<Integer, CodonTable> allCodonTables;
    private CodonTable currentCodonTable;

    private static AminoAcidManager instance;

    private AminoAcidManager() {
        initAANameMap();
    }

    public static AminoAcidManager getInstance() {
        if (instance == null) {
            try {
                setCodonTablesPath(codonTablesPath);
            } catch (IOException e) {
                handleExceptionLoading(e);
            } catch (JSONException e) {
                handleExceptionLoading(e);
            }

        }
        return instance;
    }

    /**
     * Reset the codon table to the default file,
     * and the current codon table to the default contained
     * in that file
     *
     * @return Instance of AminoAcidManager, for chaining
     */
    public static AminoAcidManager resetToDefaultCodonTables() {
        codonTablesPath = DEFAULT_CODON_TABLE_PATH;
        return getInstance();
    }

    private static void handleExceptionLoading(Exception e) {
        log.error(e);
        if (instance == null) {
            throw new IllegalStateException("No codon table present, and error loading " + codonTablesPath, e);
        }
    }

    public static AminoAcidManager setCodonTablesPath(String path) throws IOException, JSONException {
        AminoAcidManager newInstance = new AminoAcidManager();

        newInstance.loadCodonTables(path);
        codonTablesPath = path;
        instance = newInstance;
        return instance;
    }

    /**
     * Each codon translation table is identified by an integer id
     * These are specified in the file
     *
     * @param id "id" field in file
     * @return Whether setting the table was successful
     */
    public boolean setCodonTableById(int id) {
        if (allCodonTables.containsKey(id)) {
            currentCodonTable = allCodonTables.get(id);
            return true;
        } else {
            return false;
        }
    }

    /**
     * @param codon 3-letter nucleotide sequence
     * @return The amino acid represented by this codon, as
     *         decoded from the current codon table
     */
    public AminoAcid getAminoAcid(String codon) {
        return currentCodonTable.getAminoAcid(codon);
    }

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
     * @return AminoAcidSequence, or null if seqBytes == null
     */
    public AminoAcidSequence getAminoAcidSequence(byte[] seqBytes, int startPosition, Strand strand) {
        if (seqBytes == null) {
            return null;
        } else {
            String nucSequence = new String(seqBytes);
            List<AminoAcid> acids = getAminoAcids(nucSequence, strand);

            // Set the start position of this amino acid.
            return new AminoAcidSequence(strand, startPosition, acids);
        }
    }

    /**
     * Return an amino acid sequence for the input sequence.
     *
     * @param sequence
     * @param direction
     * @return
     */
    public List<AminoAcid> getAminoAcids(String sequence, Strand direction) {

        // Sequence must be divisible by 3. It is the responsibility of the
        // calling program to send a sequence properly aligned.  
        int readLength = sequence.length() / 3;
        List<AminoAcid> acids = new ArrayList<AminoAcid>(readLength);

        for (int i = 0; i <= sequence.length() - 3; i += 3) {
            String codon = sequence.substring(i, i + 3).toUpperCase();
            if (direction == Strand.NEGATIVE) {
                codon = getNucleotideComplement(codon);
            }
            AminoAcid aa = currentCodonTable.getAminoAcid(codon);
            acids.add(aa);
        }
        return acids;
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

    public Set<String> getMappingSNPs(String codon, AminoAcid mutAA) {
        Set<String> mapSNPs = new HashSet<String>();
        Set<String> SNPs = getAllSNPs(codon);
        for (String modCodon : SNPs) {
            //We use short name because all 3 stop codon have different long names,
            //and we don't care about the difference here.
            if (currentCodonTable.getAminoAcid(modCodon).equalsByName(mutAA.getShortName())) {
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
     * Load codon tables from the specified path. If any exceptions occur
     * while loading, no changes are made to this instance;
     * <p/>
     * The currentCodonTable is set to be the first codonTable in
     * the JSONArray
     *
     * @param codonTablesPath
     * @return
     */
    private synchronized void loadCodonTables(String codonTablesPath) throws IOException, JSONException {
        LinkedHashMap<Integer, CodonTable> newCodonTables = new LinkedHashMap<Integer, CodonTable>(20);

        InputStream is = AminoAcidManager.class.getResourceAsStream(codonTablesPath);
        if (is == null) {
            is = ParsingUtils.openInputStream(codonTablesPath);
        }

        JSONObject allData = readFromStream(is);
        JSONArray codonArray = allData.getJSONArray("Genetic-code-table");
        if (codonArray.length() == 0) {
            throw new JSONException("JSON File has empty array for Genetic-code-table");
        }
        CodonTable firstCodonTable = null;
        for (int ca = 0; ca < codonArray.length(); ca++) {
            CodonTable curTable = CodonTable.createFromJSON(codonArray.getJSONObject(ca));
            newCodonTables.put(curTable.getId(), curTable);
            if (firstCodonTable == null) firstCodonTable = curTable;
        }

        allCodonTables = newCodonTables;
        currentCodonTable = firstCodonTable;
    }

    private JSONObject readFromStream(InputStream is) throws JSONException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        JSONTokener tokener = new JSONTokener(reader);
        return new JSONObject(tokener);
    }

    /**
     * Initialize table of amino acid names, for easy lookup of
     * AminoAcid by symbols. This method is idempotent, only called once
     * to read name file.
     */
    private synchronized static void initAANameMap() {
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

    /**
     * Store information about current codon translation table.
     * Intended to be loaded from external resource, and then never modified.
     * To that end, collections contained here are set to be unmodifiable
     */
    public static class CodonTable {

        private final int id;
        private final List<String> names;

        private final Set<AminoAcid> starts;
        private final Map<String, AminoAcid> codonMap;

        /**
         * Get the amino acid represented by this codon
         *
         * @param codon
         * @return
         */
        public AminoAcid getAminoAcid(String codon) {
            if (codon.length() != 3) {
                throw new IllegalArgumentException("Codon must be length 3: " + codon);
            }

            AminoAcid aa = codonMap.get(codon);
            if (aa == null) {
                return AminoAcid.NULL_AMINO_ACID;
            }
            return aa;
        }

        private CodonTable(int id, List<String> names, Set<AminoAcid> starts, Map<String, AminoAcid> codonMap) {
            this.id = id;
            this.names = Collections.unmodifiableList(names);
            this.starts = Collections.unmodifiableSet(starts);
            this.codonMap = Collections.unmodifiableMap(codonMap);
        }

        private static CodonTable createFromJSON(JSONObject jsonObject) throws JSONException {
            int id = jsonObject.getInt("id");

            JSONArray jsonnames = jsonObject.getJSONArray("name");

            List<String> names = new ArrayList<String>(jsonnames.length());

            for (int nn = 0; nn < jsonnames.length(); nn++) {
                names.add(jsonnames.getString(nn));
            }

            //Data is written as several long strings which line up
            String aas = jsonObject.getString("ncbieaa");
            String startString = jsonObject.getString("sncbieaa");

            String base1 = jsonObject.getString("Base1");
            String base2 = jsonObject.getString("Base2");
            String base3 = jsonObject.getString("Base3");

            checkLengths(aas, startString, base1, base2, base3);

            Map<String, AminoAcid> codonMap = new HashMap<String, AminoAcid>(aas.length());
            Set<AminoAcid> starts = new HashSet<AminoAcid>(aas.length());

            for (int cc = 0; cc < aas.length(); cc++) {
                String codon = base1.substring(cc, cc + 1) + base2.substring(cc, cc + 1) + base3.substring(cc, cc + 1);
                AminoAcid aa = AANameMap.get(aas.substring(cc, cc + 1));

                codonMap.put(codon, aa);

                if (startString.charAt(cc) == 'M') {
                    starts.add(aa);
                }
            }

            return new CodonTable(id, names, starts, codonMap);

        }

        private static void checkLengths(String... values) throws JSONException {
            int length = values[0].length();
            for (int v = 1; v < values.length; v++) {
                if (values[v].length() != length) {
                    String msg = "Amino acid and codon strings must all be the same length.";
                    msg += "Expected length " + length + ", found length " + values[v].length();
                    throw new JSONException(msg);
                }
            }
        }

        public int getId() {
            return id;
        }

        public String getDisplayName() {
            return names.get(0);
        }

        public Set<AminoAcid> getStarts() {
            return starts;
        }
    }

}
