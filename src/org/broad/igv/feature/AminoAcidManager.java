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
package org.broad.igv.feature;

import com.google.common.base.Objects;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import com.google.gson.*;
import org.apache.log4j.Logger;
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

    private static final String[] BASE_SEQUENCES = {"TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
            "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
            "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"};


    static final String DEFAULT_CODON_TABLE_PATH = "resources/geneticCode.json";

    static final String DEFAULT_TRANS_TABLE_PATH = "resources/defaultTranslationTables.json";

    //ID of the "standard" translation table
    public static final int STANDARD_TABLE_ID = 1;

    private static final String DEFAULT_CHROMO_KEY = "default";

    private LinkedHashMap<CodonTableKey, CodonTable> allCodonTables = new LinkedHashMap<CodonTableKey, CodonTable>(20);
    private CodonTable currentCodonTable;

    private static Table<String, String, CodonTableKey> genomeChromoTable = TreeBasedTable.create();

    private static AminoAcidManager instance;

    private AminoAcidManager() {
        initAANameMap();
        try {
            loadDefaultTranslationTables();
        } catch (JsonParseException e) {
            log.error(e);
        }
    }

    public static AminoAcidManager getInstance() {
        if (instance == null) {
            try {
                AminoAcidManager newInstance = new AminoAcidManager();
                newInstance.loadCodonTables(DEFAULT_CODON_TABLE_PATH);
                instance = newInstance;
            } catch (IOException e) {
                handleExceptionLoading(e);
            } catch (JsonParseException e) {
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
        instance = null;
        return getInstance();
    }

    private static void handleExceptionLoading(Exception e) {
        log.error(e);
        if (instance == null) {
            throw new IllegalStateException("No codon table present, and error loading " + DEFAULT_CODON_TABLE_PATH, e);
        }
    }

    /**
     * Removes all codon tables.
     * Mainly for testing
     */
    synchronized void clear() {
        allCodonTables.clear();
        currentCodonTable = null;
    }

    /**
     * Each codon translation table is identified by an integer id
     * These are specified in the file. We specify a table
     * by filename/id combination
     *
     * @param codonTablePath
     * @param id
     * @return Whether setting the table was successful
     */
    public boolean setCodonTable(String codonTablePath, int id) {
        CodonTableKey key = new CodonTableKey(codonTablePath, id);
        return setCodonTable(key);
    }

    public boolean setCodonTable(CodonTableKey key) {
        if (allCodonTables.containsKey(key)) {
            currentCodonTable = allCodonTables.get(key);
            return true;
        } else {
            return false;
        }
    }

    /**
     * @param codon 3-letter nucleotide sequence
     * @return The amino acid represented by this codon, as
     * decoded from the current codon table
     */
    public AminoAcid getAminoAcid(String codon) {
        return currentCodonTable.getAminoAcid(codon);
    }


    /**
     * Return a list of amino acids for the input sequence of nucleotides
     *
     * @param direction
     * @param sequence
     * @return
     */
    List<AminoAcid> getAminoAcids(Strand direction, String sequence) {

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
     * Get the amino acid sequence for an interval.
     * Assumptions and conventions
     * <p/>
     * The start and end positions are on the positive strand
     * irrespective of the read direction.
     * <p/>
     * Reading will begin from the startPosition if strand == POSITIVE, endPosition if NEGATIVE
     *
     * @param strand
     * @param startPosition
     * @param seqBytes
     * @return AminoAcidSequence, or null if seqBytes == null
     */
    public synchronized AminoAcidSequence getAminoAcidSequence(Strand strand, int startPosition, byte[] seqBytes) {
        if (seqBytes == null) {
            return null;
        } else {
            String nucSequence = new String(seqBytes);
            List<AminoAcid> acids = getAminoAcids(strand, nucSequence);
            return new AminoAcidSequence(strand, startPosition, acids, currentCodonTable.getKey());
        }
    }

    /**
     * Return a list of amino acids for the input nucleotides
     *
     * @param strand
     * @param startPosition
     * @param nucleotides
     * @return
     */
    public synchronized AminoAcidSequence getAminoAcidSequence(Strand strand, int startPosition, String nucleotides) {
        if (nucleotides == null) {
            return null;
        } else {
            List<AminoAcid> acids = getAminoAcids(strand, nucleotides);
            return new AminoAcidSequence(strand, startPosition, acids, currentCodonTable.getKey());
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
     * while loading, no changes are made to this instance.
     * <p/>
     * Note that the new codon tables are ADDED to the existing tables
     * <p/>
     * The currentCodonTable is set to be the codonTable with id = defaultid if present
     * If not, the first one in the array is set as default
     *
     * @param codonTablesPath
     * @return
     */
    synchronized void loadCodonTables(String codonTablesPath) throws IOException, JsonParseException {
        LinkedHashMap<CodonTableKey, CodonTable> newCodonTables = new LinkedHashMap<CodonTableKey, CodonTable>(20);
        CodonTable defaultCodonTable = null;

        InputStream is = AminoAcidManager.class.getResourceAsStream(codonTablesPath);
        if (is == null) {
            is = ParsingUtils.openInputStream(codonTablesPath);
        }

        if (codonTablesPath.endsWith(".json")) {
            JsonObject allData = readJSONFromStream(is);
            int defaultId = -1;
            defaultId = allData.get("defaultid").getAsInt();
            JsonArray codonArray = allData.get("Genetic-code-table").getAsJsonArray();
            if (codonArray.size() == 0) {
                throw new JsonParseException("JSON File has empty array for Genetic-code-table");
            }
            for (int ca = 0; ca < codonArray.size(); ca++) {
                CodonTable curTable = CodonTable.createFromJSON(codonTablesPath, codonArray.get(ca).getAsJsonObject());
                newCodonTables.put(curTable.getKey(), curTable);
                if (defaultCodonTable == null || curTable.getId() == defaultId) {
                    defaultCodonTable = curTable;
                }
            }
        } else {
            throw new IllegalArgumentException("Unknown file type, must be .json");
        }

        allCodonTables.putAll(newCodonTables);
        currentCodonTable = defaultCodonTable;
    }

//    private static JsonObject readJSONFromStream(InputStream is) throws JsonParseException {
//        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
//        JSONTokener tokener = new JSONTokener(reader);
//        return new JsonObject(tokener);
//    }

    private static JsonObject readJSONFromStream(InputStream is) {
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        JsonParser parser = new JsonParser();
        return parser.parse(reader).getAsJsonObject();
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

    public Collection<CodonTable> getAllCodonTables() {
        return Collections.unmodifiableCollection(allCodonTables.values());
    }

    public CodonTable getCodonTable() {
        return currentCodonTable;
    }

    private static void loadDefaultTranslationTables() throws JsonParseException {
        InputStream is = AminoAcidManager.class.getResourceAsStream(DEFAULT_TRANS_TABLE_PATH);
        JsonObject allData = readJSONFromStream(is);
        JsonArray organisms = allData.get("organisms").getAsJsonArray();

        for (int ind = 0; ind < organisms.size(); ind++) {
            JsonObject obj = organisms.get(ind).getAsJsonObject();

            //Process each translation table setting
            String genomeId = obj.get("genomeId").getAsString();

            String codonTablePath = DEFAULT_CODON_TABLE_PATH;
            try {
                Object tmpPath = obj.get("codonTablePath");
                if (tmpPath != null && tmpPath != JsonNull.INSTANCE && tmpPath instanceof String) {
                    codonTablePath = (String) tmpPath;
                }
            } catch (JsonParseException e) {
                log.error("No codon table path found in " + DEFAULT_TRANS_TABLE_PATH + ". Using default: " + codonTablePath);
            }

            JsonObject chromosomes = obj.get("chromosomes").getAsJsonObject();
            Iterator<Map.Entry<String, JsonElement>> iterator = chromosomes.entrySet().iterator();
            while (iterator.hasNext()) {
                Map.Entry<String, JsonElement> entry = iterator.next();
                String chromoName = entry.getKey();
                int id = entry.getValue().getAsInt();
                CodonTableKey key = new CodonTableKey(codonTablePath, id);
                genomeChromoTable.put(genomeId, chromoName, key);
            }

        }

    }

//    /**
//     * Load the default codon table for the given genome and chromosome.
//     * We check the given name, alias, and finally use the default for the specified
//     * genome.
//     *
//     * @param genome
//     * @param chrName
//     */
//    public void loadDefaultCodonTable(Genome genome, String chrName) {
//        Map<String, CodonTableKey> chrMap = genomeChromoTable.row(genome.getId());
//        String[] tryChromos = new String[]{
//                chrName, genome.getChromosomeAlias(chrName), DEFAULT_CHROMO_KEY
//        };
//        for (String tryChromo : tryChromos) {
//            if (chrMap.containsKey(tryChromo)) {
//                setCodonTable(chrMap.get(tryChromo));
//                return;
//            }
//        }
//    }

    public static class CodonTableKey {

        private final String sourcePath;
        private final int id;

        private CodonTableKey(String sourcePath, int id) {
            this.sourcePath = sourcePath;
            this.id = id;
        }

        @Override
        public boolean equals(Object object) {
            if (object instanceof CodonTableKey) {
                CodonTableKey other = (CodonTableKey) object;
                return this.id == other.id &&
                        Objects.equal(this.sourcePath, other.sourcePath);
            }
            return false;
        }

        @Override
        public int hashCode() {
            return Objects.hashCode(this.sourcePath, this.id);
        }

        public int getId() {
            return id;
        }
    }

    /**
     * Store information about current codon translation table.
     * Intended to be loaded from external resource, and then never modified.
     * To that end, collections contained here are set to be unmodifiable
     */
    public static class CodonTable {

        private final CodonTableKey key;
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

        private CodonTable(String path, int id, List<String> names, Set<AminoAcid> starts, Map<String, AminoAcid> codonMap) {
            this.key = new CodonTableKey(path, id);
            this.names = Collections.unmodifiableList(names);
            this.starts = Collections.unmodifiableSet(starts);
            this.codonMap = Collections.unmodifiableMap(codonMap);
        }

        private static CodonTable createFromJSON(String sourcePath, JsonObject jsonObject) throws JsonParseException {
            int id = jsonObject.get("id").getAsInt();

            JsonArray jsonnames = jsonObject.get("name").getAsJsonArray();
            List<String> names = new ArrayList<String>(jsonnames.size());
            for (int nn = 0; nn < jsonnames.size(); nn++) {
                names.add(jsonnames.get(nn).getAsString());
            }

            //Data is written as several long strings which line up
            String aas = jsonObject.get("ncbieaa").getAsString();
            String startString = jsonObject.get("sncbieaa").getAsString();

            return build(sourcePath, id, names, aas, startString);
        }

        private static CodonTable build(String sourcePath, int id, List<String> names, String aas, String startString) {

            String base1 = BASE_SEQUENCES[0];
            String base2 = BASE_SEQUENCES[1];
            String base3 = BASE_SEQUENCES[2];

            checkLengths(base1, base2, base3, aas, startString);


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

            return new CodonTable(sourcePath, id, names, starts, codonMap);
        }

        private static void checkLengths(String... values) {
            int length = values[0].length();
            assert length == 64;
            for (int v = 1; v < values.length; v++) {
                if (values[v].length() != length) {
                    String msg = "Amino acid and codon strings must all be the same length.";
                    msg += "Expected length " + length + ", found length " + values[v].length();
                    throw new InputMismatchException(msg);
                }
            }
        }

        public int getId() {
            return key.id;
        }

        public String getDisplayName() {
            return names.get(0);
        }

        public Set<AminoAcid> getStarts() {
            return starts;
        }

        Map<String, AminoAcid> getCodonMap() {
            return codonMap;
        }

        @Override
        public boolean equals(Object object) {
            if (object instanceof CodonTable) {
                CodonTable other = (CodonTable) object;
                return Objects.equal(this.key, other.key) &&
                        Objects.equal(this.names, other.names) &&
                        Objects.equal(this.starts, other.starts) &&
                        Objects.equal(this.codonMap, other.codonMap);
            }
            return false;
        }

        @Override
        public int hashCode() {
            return Objects.hashCode(this.key.id, this.key.sourcePath, this.names, this.starts, this.codonMap);
        }

        public CodonTableKey getKey() {
            return key;
        }
    }

}
