package org.broad.igv.feature.genome;

import org.broad.igv.Globals;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.genome.GenomeDescriptor;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.net.URL;
import java.util.*;

import static org.broad.igv.prefs.Constants.BACKUP_GENOMES_SERVER_URL;
import static org.broad.igv.prefs.Constants.GENOMES_SERVER_URL;

/**
 * Singleton class for managing the list of genomes presented in the command bar combo box.  Also has methods for
 * searching the hosted genome list.
 */
public class HostedGenomes {

    private static Logger log = LogManager.getLogger(HostedGenomes.class);

    private static String[] legacyColumns = {
            "common name",
            "url",
            "assembly"
    };


    private static CopyOnWriteArrayList<GenomeDescriptor> records;

    private static Map<String, GenomeDescriptor> hostedGenomesMap = null;

    public static List<GenomeDescriptor> getRecords() {
        if (records == null) {
            records = new CopyOnWriteArrayList<>(readRecords());
        }
        return records;
    }


    public static GenomeDescriptor getGenomeTableRecord(String genomeId) {
        if (hostedGenomesMap == null) {
            hostedGenomesMap = new HashMap<>();
            for (GenomeDescriptor record : getRecords()) {
                hostedGenomesMap.put(record.getId(), record);
            }
        }
        return hostedGenomesMap.get(genomeId);
    }

private static List<GenomeDescriptor> readRecords() {

    records = new ArrayList<>();

    final IGVPreferences preferences = PreferencesManager.getPreferences();
    final String genomesServerURL = preferences.get(GENOMES_SERVER_URL);
    final String backupGenomesServerURL = preferences.get(BACKUP_GENOMES_SERVER_URL);
    final String genarkURL = "https://hgdownload.soe.ucsc.edu/hubs/UCSC_GI.assemblyHubList.txt";

    List<String> errors = new ArrayList<>();

    // IGV hosted genome list
    boolean genomeListLoaded = loadGenomeList(genomesServerURL, "assembly", errors);
    if (!genomeListLoaded) {
        log.error("Error loading genome list from: " + genomesServerURL);
        errors.add("Error loading genome list from: " + genomesServerURL);
        // Try backup server
        if (!loadGenomeList(backupGenomesServerURL, "assembly", errors)) {
            errors.add("Error loading genome list from: " + backupGenomesServerURL);
        }
    }

    // UCSC Genark hosted genome list
    if (!loadGenomeList(genarkURL, "assembly", errors)) {
        log.error("Error connecting to UCSC Genark server URL: " + genarkURL);
        errors.add("Error connecting to UCSC Genark server: " + genarkURL);
    }

    if (!errors.isEmpty()) {
        StringBuilder sb = new StringBuilder();
        for (String error : errors) {
            sb.append(error).append("\n");
        }
        MessageUtils.showMessage(sb.toString());
    }

    return records;
}

private static boolean loadGenomeList(String url, String idColumn, List<String> errors) {
    try {
        String genomeListContent = HttpUtils.getInstance().getContentsAsString(new URL(url));
        List<String> genomeListLines = Arrays.asList(genomeListContent.split("\\r?\\n"));
        String[] headers = parseHeaders(genomeListLines);
        parseRecords(genomeListLines, headers, idColumn);
        return true;
    } catch (Exception e) {
        log.error("Error loading genome list from: " + url, e);
        errors.add("Error loading genome list from: " + url + "   (" + e.getMessage() + ")");
        return false;
    }
}

    private static String[] parseHeaders(List<String> genomeListLines) {

        String [] headers = null;

        // Find last line starting with "#"
        String lastHeaderLine = null;
        for (int i = 0; i < genomeListLines.size(); i++) {
            String line = genomeListLines.get(i).trim();
            if (line.startsWith("#")) {
                lastHeaderLine = line;
            } else {
                break;
            }
        }

        if (lastHeaderLine != null) {
            String[] tokens = Arrays.stream(Globals.tabPattern.split(lastHeaderLine.substring(1))).map(h -> h.trim()).toArray(String[]::new);
            if (tokens.length >= 3) {
                headers = tokens;
            } else {
                headers = legacyColumns;
            }
        } else {
            headers = legacyColumns;
        }

        return headers;
    }


    private static void parseRecords(List<String> genomeListLines, String [] headers, String idColumn) {

        for (String line : genomeListLines) {

            if (line.startsWith("<Server-Side>") || line.startsWith("#")) {
                // Skip header lines and server-side comments
            } else {
                line = line.trim();
                String[] values = Globals.tabPattern.split(line.trim());

                Map<String, String> attributes = new HashMap<>();
                for (int i = 0; i < headers.length; i++) {
                    attributes.put(headers[i], values[i]);
                }
                String id = attributes.get(idColumn);
                String displayableName = attributes.get("common name");
                String path = attributes.get("url");
                records.add(new GenomeDescriptor(displayableName, path, id, attributes));
            }
        }
    }
}
