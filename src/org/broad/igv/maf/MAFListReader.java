package org.broad.igv.maf;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Implementation of MAFReader for MAF files that are split by chromosome (1 file per chromosome).  Requires
 * a 2 column mapping file,  chr -> maf file.
 *
 * @author Jim Robinson
 * @date 4/16/12
 */
public class MAFListReader implements MAFReader {

    private static Logger log = Logger.getLogger(MAFListReader.class);

    List<String> chrNames;

    // Map of chr name -> MAF file path
    Map<String, String> filenameMap;

    // Map of chr name -> MAFLocalReader
    Map<String, MAFLocalReader> readerMap;

    public MAFListReader(String mappingFile) throws IOException {

        loadDictionaryFile(mappingFile);
        readerMap = new HashMap<String, MAFLocalReader>();

    }

    private void loadDictionaryFile(String mappingFile) throws IOException {

        chrNames = new ArrayList<String>();
        // Map of chr name -> MAF file path
        filenameMap = new HashMap();


        // TODO -- this is a very common pattern, a 2 column file representing a dictionary.
        // TODO    Create some utility method for this
        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(mappingFile);
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("#")) continue;
                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                if (tokens.length != 2) {
                    log.info("Skipping line: " + nextLine);
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

    public MAFTile loadTile(String chr, int start, int end, List<String> species) {

        MAFLocalReader reader = getReader(chr);
        return reader == null ? null : reader.loadTile(chr, start, end, species);
    }

    private MAFLocalReader getReader(final String chr) {
        MAFLocalReader reader = readerMap.get(chr);
        if (reader == null) {
            final String path = filenameMap.get(chr);
            if (path == null) {
                log.info("No MAF file found for chromosome: " + chr);
            } else {
                Runnable runnable = new Runnable() {
                    public void run() {
                        try {
                            MAFLocalReader reader = new MAFLocalReader(path);
                            readerMap.put(chr, reader);
                            IGV.getInstance().repaintDataAndHeaderPanels();
                        } catch (Exception e) {
                            log.error("Error loading MAF reader (" + path + "):  ", e);
                            MessageUtils.showMessage("Error loading MAF file: " + e.getMessage());
                        }
                    }
                };
                LongRunningTask.submit(runnable);
            }

        }
        return reader;
    }


    public List<String> getChrNames() {
        return chrNames;
    }
}
