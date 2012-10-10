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

import com.google.common.base.Objects;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import org.apache.log4j.Logger;
import org.bouncycastle.asn1.*;
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
        } catch (JSONException e) {
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
     *         decoded from the current codon table
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
    synchronized void loadCodonTables(String codonTablesPath) throws IOException, JSONException {
        LinkedHashMap<CodonTableKey, CodonTable> newCodonTables = new LinkedHashMap<CodonTableKey, CodonTable>(20);
        CodonTable defaultCodonTable = null;

        InputStream is = AminoAcidManager.class.getResourceAsStream(codonTablesPath);
        if (is == null) {
            is = ParsingUtils.openInputStream(codonTablesPath);
        }

        if (codonTablesPath.endsWith(".json")) {
            JSONObject allData = readJSONFromStream(is);
            int defaultId = -1;
            try {
                defaultId = allData.getInt("defaultid");
            } catch (JSONException e) {
                //pass;
            }
            JSONArray codonArray = allData.getJSONArray("Genetic-code-table");
            if (codonArray.length() == 0) {
                throw new JSONException("JSON File has empty array for Genetic-code-table");
            }
            for (int ca = 0; ca < codonArray.length(); ca++) {
                CodonTable curTable = CodonTable.createFromJSON(codonTablesPath, codonArray.getJSONObject(ca));
                newCodonTables.put(curTable.getKey(), curTable);
                if (defaultCodonTable == null || curTable.getId() == defaultId) {
                    defaultCodonTable = curTable;
                }
            }
        } else if (codonTablesPath.endsWith(".asn1") || codonTablesPath.endsWith(".val")) {
            ASN1InputStream ASNis = new ASN1InputStream(is);
            ASN1Primitive obj = ASNis.readObject();
            ASN1Set set = (ASN1Set) obj;
            //Array of different genetic code tables
            ASN1Encodable[] codonArray = set.toArray();
            if (codonArray.length == 0) {
                throw new RuntimeException("ASN1 File has empty array for Genetic-code-table");
            }
            for (ASN1Encodable aCodonArray : codonArray) {
                CodonTable curTable = CodonTable.createFromASN1(codonTablesPath, aCodonArray);
                newCodonTables.put(curTable.getKey(), curTable);
                if (defaultCodonTable == null) {
                    defaultCodonTable = curTable;
                }
            }
        } else {
            throw new IllegalArgumentException("Unknown file type, must be .json or .asn1");
        }

        allCodonTables.putAll(newCodonTables);
        currentCodonTable = defaultCodonTable;
    }

    private static JSONObject readJSONFromStream(InputStream is) throws JSONException {
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

    public Collection<CodonTable> getAllCodonTables() {
        return Collections.unmodifiableCollection(allCodonTables.values());
    }

    public CodonTable getCodonTable() {
        return currentCodonTable;
    }

    private static void loadDefaultTranslationTables() throws JSONException {
        InputStream is = AminoAcidManager.class.getResourceAsStream(DEFAULT_TRANS_TABLE_PATH);
        JSONObject allData = readJSONFromStream(is);
        JSONArray organisms = allData.getJSONArray("organisms");

        for (int ind = 0; ind < organisms.length(); ind++) {
            JSONObject obj = organisms.getJSONObject(ind);

            //Process each translation table setting
            String genomeId = obj.getString("genomeId");

            String codonTablePath = DEFAULT_CODON_TABLE_PATH;
            try {
                Object tmpPath = obj.get("codonTablePath");
                if (tmpPath != null && tmpPath != JSONObject.NULL && tmpPath instanceof String) {
                    codonTablePath = (String) tmpPath;
                }
            } catch (JSONException e) {
                log.error("No codon table path found in " + DEFAULT_TRANS_TABLE_PATH + ". Using default: " + codonTablePath);
            }

            JSONObject chromosomes = obj.getJSONObject("chromosomes");
            Iterator<String> iterator = chromosomes.keys();
            while (iterator.hasNext()) {
                String chromoName = iterator.next();
                int id = chromosomes.getInt(chromoName);
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

        private static CodonTable createFromJSON(String sourcePath, JSONObject jsonObject) throws JSONException {
            int id = jsonObject.getInt("id");

            JSONArray jsonnames = jsonObject.getJSONArray("name");
            List<String> names = new ArrayList<String>(jsonnames.length());
            for (int nn = 0; nn < jsonnames.length(); nn++) {
                names.add(jsonnames.getString(nn));
            }

            //Data is written as several long strings which line up
            String aas = jsonObject.getString("ncbieaa");
            String startString = jsonObject.getString("sncbieaa");

            return build(sourcePath, id, names, aas, startString);
        }

        private static CodonTable createFromASN1(String sourcePath, ASN1Encodable asn1Encodable) throws IOException {
            byte[] data = asn1Encodable.toASN1Primitive().getEncoded();
            ASN1InputStream iASNis = new ASN1InputStream(data);
            ASN1Primitive prim = iASNis.readObject();
            ASN1Set iset = (ASN1Set) prim;

            //Set of fields of each table
            ASN1TaggedObject[] taggedObjects = getTaggedObjects(iset.toArray());
            int index = 0;
            int tagNo = taggedObjects[index].getTagNo();
            List<String> names = new ArrayList<String>(2);
            while (tagNo == 0) {
                names.add(getAsString(taggedObjects[index].getObject()));
                tagNo = taggedObjects[++index].getTagNo();
            }

            int id = ((DERInteger) taggedObjects[index++].getObject()).getValue().intValue();
            String aas = getAsString(taggedObjects[index++].getObject());
            String startString = getAsString(taggedObjects[index++].getObject());

            return build(sourcePath, id, names, aas, startString);
        }

        private static String getAsString(ASN1Object object) {
            return ((ASN1String) object).getString();
        }

        private static ASN1TaggedObject[] getTaggedObjects(ASN1Encodable[] encodables) {
            ASN1TaggedObject[] taggedObjects = new ASN1TaggedObject[encodables.length];
            for (int ii = 0; ii < encodables.length; ii++) {
                taggedObjects[ii] = (ASN1TaggedObject) encodables[ii];
            }
            return taggedObjects;
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
