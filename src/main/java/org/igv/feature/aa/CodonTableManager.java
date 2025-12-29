/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.feature.aa;

import org.igv.logging.*;
import org.igv.feature.genome.GenomeManager;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * @author jrobinso
 */
public class CodonTableManager {

    private static final Logger log = LogManager.getLogger(CodonTableManager.class);

    /**
     * The genetic codes, mapping codons -> amino acides
     */
    public static final String DEFAULT_CODON_TABLE_PATH = "resources/geneticCode.json";

    /**
     * Mappings of organism -> genetic code, organism in IGV terms is defined by genomeID + chr name.  This is a bit
     * odd, but by convention the Mitochondria organelle is often represented as a "chromosome" of a genome assembly,
     * thus  hg19 / chrM can have a different genetic code than the genomic chromosomes.
     */
    static final String DEFAULT_TRANS_TABLE_PATH = "resources/defaultTranslationTables.json";

    //ID of the "standard" translation table
    public static final int STANDARD_TABLE_ID = 1;

    private Map<String, CodonTableMap> genomeChromoTable = new HashMap<>();

    private LinkedHashMap<Integer, CodonTable> allCodonTables = new LinkedHashMap<>(20);

    private CodonTable defaultCodonTable;

    /**
     * Explicitly set codon table (e.g. by user action), overrides default.
     */
    private CodonTable currentCodonTable;

    private static CodonTableManager instance;

    public static synchronized CodonTableManager getInstance() {
        if (instance == null) {
            AminoAcidManager.initAANameMap();    // TODO -- where should this be done?
            instance = new CodonTableManager();
        }
        return instance;
    }

    private CodonTableManager() {
        init();
    }

    private void init() {
        try {
            loadDefaultTranslationTables();
            loadCodonTables(DEFAULT_CODON_TABLE_PATH);
        } catch (IOException e) {
            handleExceptionLoading(e);
        } catch (Exception e) {
            log.error(e);
        }
    }


    private static void handleExceptionLoading(Exception e) {
        log.error(e);
        if (instance == null) {
            throw new IllegalStateException("No codon table present, and error loading " + DEFAULT_CODON_TABLE_PATH, e);
        }
    }

    /**
     * Reset the codon table to the default file, and the current codon table to the default contained
     * in that file.  Useful for unit testing.
     *
     * @return Instance of AminoAcidManager, for chaining
     */
    public void resetToDefaults() {
        currentCodonTable = null;
        genomeChromoTable = new HashMap<>();
        allCodonTables = new LinkedHashMap<>(20);
        init();
    }

    /**
     * Return the codon table for the given organelle, defined by a genomeID and chr name.  This accomodates the
     * common, if biologically suspect, convention of treating the mitochondria organelle as a "chromosome" of an organism
     *
     * @param genomeID
     * @param chr
     * @return
     */
    public CodonTable getCodonTableForChromosome(String genomeID, String chr) {

        if (currentCodonTable != null) {
            return currentCodonTable;
        } else {
            CodonTableMap map = genomeChromoTable.get(genomeID);
            if (map == null) {
                return defaultCodonTable;
            } else {
                Integer tableID = map.getTableIdForChr(chr);
                return allCodonTables.get(tableID);
            }
        }
    }

    /**
     * Convenience method
     *
     * @param chr
     * @return
     */
    public CodonTable getCodonTableForChromosome(String chr) {
        return getCodonTableForChromosome(GenomeManager.getInstance().getGenomeId(), chr);
    }

    /**
     * Return all codon tables.  Used to popuplate menus.
     *
     * @return
     */
    public Collection<CodonTable> getAllCodonTables() {
        return Collections.unmodifiableCollection(allCodonTables.values());
    }

    /**
     * Explicitly set the codon table (e.g. by user action), overrides default.
     */
    public void setCurrentCodonTable(CodonTable codonTable) {
        currentCodonTable = codonTable;
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
    synchronized void loadCodonTables(String codonTablesPath) throws IOException {

        LinkedHashMap<Integer, CodonTable> newCodonTables = new LinkedHashMap<>(20);

        InputStream is = AminoAcidManager.class.getResourceAsStream(codonTablesPath);
        if (is == null) {
            is = ParsingUtils.openInputStream(codonTablesPath);
        }

        if (codonTablesPath.endsWith(".json")) {
            StringBuilder sb = new StringBuilder();
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(is))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    sb.append(line);
                }
            }
            org.json.JSONObject allData = new org.json.JSONObject(sb.toString());
            int defaultId = allData.getInt("defaultid");
            org.json.JSONArray codonArray = allData.getJSONArray("Genetic-code-table");
            if (codonArray.length() == 0) {
                throw new RuntimeException("JSON File has empty array for Genetic-code-table");
            }
            for (int ca = 0; ca < codonArray.length(); ca++) {
                org.json.JSONObject codonObj = codonArray.getJSONObject(ca);
                CodonTable curTable = CodonTable.createFromJSON(codonTablesPath, codonObj);
                newCodonTables.put(curTable.getId(), curTable);
                if (defaultCodonTable == null || curTable.getId() == defaultId) {
                    defaultCodonTable = curTable;
                }
            }
        } else {
            throw new IllegalArgumentException("Unknown file type, must be .json");
        }

        allCodonTables.putAll(newCodonTables);
        //currentCodonTable = defaultCodonTable;

        is.close();
    }

    private void loadDefaultTranslationTables() throws org.json.JSONException {
        InputStream is = CodonTableManager.class.getResourceAsStream(DEFAULT_TRANS_TABLE_PATH);
        org.json.JSONObject allData = new org.json.JSONObject(new java.util.Scanner(is).useDelimiter("\\A").next());
        org.json.JSONArray organisms = allData.getJSONArray("organisms");
        for (int ind = 0; ind < organisms.length(); ind++) {
            org.json.JSONObject obj = organisms.getJSONObject(ind);
            String genomeId = obj.getString("genomeId");
            genomeChromoTable.put(genomeId, new CodonTableMap(obj));
        }
    }

    public CodonTable getDefaultCodonTable() {
        return defaultCodonTable;
    }

    public CodonTable getCurrentCodonTable() {
        return currentCodonTable;
    }

    public CodonTable getCodonTableByID(Integer id) {
        return allCodonTables.get(id);
    }


    /**
     * Maps chromosome names to codon tables for a specific genome.   The main purpose of this is to allow the
     * treatment of the mitochondria sequences as a "chromosome" for a genome, biologically wrong but convenient
     * for genome browsers.
     */
    static class CodonTableMap {

        String genomeID;
        Integer defaultID;
        Map<String, Integer> chromosomeIDs;

        CodonTableMap(org.json.JSONObject obj) {
            genomeID = obj.getString("genomeId");
            chromosomeIDs = new HashMap<>();

            org.json.JSONObject chromosomes = obj.getJSONObject("chromosomes");
            defaultID = chromosomes.getInt("default");
            Iterator<String> keys = chromosomes.keys();
            while (keys.hasNext()) {
                String chromoName = keys.next();
                if (!chromoName.equals("default")) {
                    int id = chromosomes.getInt(chromoName);
                    chromosomeIDs.put(chromoName, id);
                }
            }
        }

        Integer getTableIdForChr(String chr) {
            if (chromosomeIDs.containsKey(chr)) {
                return chromosomeIDs.get(chr);
            } else {
                return defaultID;
            }
        }

    }

}
