package org.igv.maf;

import org.igv.logging.*;
import org.igv.Globals;
import org.igv.feature.genome.Genome;
import org.igv.ui.IGV;
import org.igv.ui.util.MessageUtils;
import org.igv.util.FileUtils;
import org.igv.util.LongRunningTask;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Implementation of MAFReader for MAF files that are split by chromosome (1 file per chromosome).  Requires
 * a 2 column mapping file,  chr -> maf file.
 *
 * @author Jim Robinson
 * @date 4/16/12
 */
public class MAFListReader implements MAFReader {

    private static Logger log = LogManager.getLogger(MAFListReader.class);

    List<String> chrNames;

    private String refId;

    /**
     * Species (sequences) represented in this maf file.
     */
    private List<String> species;

    /**
     * A lookup table for species (id -> name).  This is optional
     */
    private Map<String, String> speciesNames;


    // Map of chr name -> MAF file path
    Map<String, String> filenameMap;

    // Map of chr name -> MAFLocalReader
    Map<String, MAFParser> readerMap;

    public MAFListReader(String mappingFile) throws IOException {

        loadDictionaryFile(mappingFile);
        loadSpeciesNames(mappingFile);
        readerMap = new HashMap<String, MAFParser>();

    }

    public String getRefId() {
        return refId;
    }

    private void loadDictionaryFile(String mappingFile) throws IOException {

        chrNames = new ArrayList<String>();
        // Map of chr name -> MAF file path
        filenameMap = new HashMap();

        BufferedReader br = null;
        try {
            br = ParsingUtils.openBufferedReader(mappingFile);
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("#")) continue;
                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                if (tokens.length != 2) {
                    log.warn("Skipping line: " + nextLine);
                } else {
                    String chr = tokens[0];
                    String fname = tokens[1];
                    String fullPath = FileUtils.getAbsolutePath(fname, mappingFile);
                    filenameMap.put(chr, fullPath);
                    chrNames.add(chr);
                }
            }
        } finally {
            if (br != null) br.close();
        }
    }

    /**
     * Load a file containing species IDs and names.
     *
     * @param path
     */
    private void loadSpeciesNames(String path) {

        InputStream is = null;
        species = new ArrayList<String>();
        speciesNames = new LinkedHashMap<String, String>();

        try {
            String speciesPath = path + ".species";
            if (FileUtils.resourceExists(speciesPath)) {
                is = ParsingUtils.openInputStream(speciesPath);
            } else {
                // No file, use default
                is = MAFUtils.class.getResourceAsStream("species.properties");
            }

            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("#ref")) {
                    String[] tokens = Globals.equalPattern.split(nextLine);
                    refId = tokens[1];
                } else {
                    String[] tokens = Globals.equalPattern.split(nextLine);
                    if (tokens.length == 2) {
                        String id = tokens[0];
                        String name = tokens[1];
                        species.add(id);
                        speciesNames.put(id, name);
                    } else {
                        //log.warn("Skipping line: " + nextLine);
                    }
                }

            }

        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {

                }
            }
        }
    }

//
//    public MAFTile loadTile(String chr, int start, int end, List<String> species) {
//
//        MAFLocalReader reader = getReader(chr);
//        return reader == null ? null : reader.loadTile(chr, start, end, species);
//    }

    @Override
    public List<MultipleAlignmentBlock> loadAlignments(String chr, int start, int end) throws IOException {
        MAFReader reader = getReader(chr);
        return reader == null ? null : reader.loadAlignments(chr, start, end);
    }

    private MAFParser getReader(final String chr) {
        MAFParser reader = readerMap.get(chr);
        if (reader == null) {
            final String path = filenameMap.get(chr);
            if (path == null) {
                log.info("No MAF file found for chromosome: " + chr);
            } else {
                try {
                    reader = new MAFParser(path);
                    readerMap.put(chr, reader);
                } catch (Exception e) {
                    log.error("Error loading MAF reader (" + path + "):  ", e);
                    MessageUtils.showMessage("Error loading MAF file: " + e.getMessage());
                }
            }
        }
        return reader;
    }


    public List<String> getChrNames() {
        return chrNames;
    }

    @Override
    public String getSpeciesName(String speciesId) {
        if (speciesNames != null && speciesNames.containsKey(speciesId)) {
            return speciesNames.get(speciesId);
        } else {
            return speciesId;
        }
    }

    public Collection<String> getSpecies() {
        return species;
    }
}
