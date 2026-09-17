package org.igv.feature.genome;

import org.igv.Globals;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.genome.GenomeListItem;
import org.igv.ui.util.MessageUtils;
import org.igv.util.HttpUtils;

import java.net.URL;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CopyOnWriteArrayList;

import static org.igv.prefs.Constants.BACKUP_GENOMES_SERVER_URL;
import static org.igv.prefs.Constants.GENOMES_SERVER_URL;

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


    private static CopyOnWriteArrayList<GenomeListItem> records;

    private static Map<String, GenomeListItem> hostedGenomesMap = null;

    /**
     * IDs of genomes listed by the IGV genome server (the GENOMES_SERVER_URL preference, or its backup).  These
     * genomes can be restored from their ID alone, so sessions reference them by ID rather than by an expanded
     * genome definition.
     */
    private static Set<String> igvHostedIds = ConcurrentHashMap.newKeySet();

    private static boolean igvHostedIdsLoaded = false;

    public static synchronized List<GenomeListItem> getRecords() {
        if (records == null) {
            records = new CopyOnWriteArrayList<>(readRecords());
        }
        return records;
    }


    /**
     * Return true if the genome ID is that of a genome hosted by the IGV genome server.  Genomes from other sources,
     * including the UCSC GenArk list, return false.
     * <p>
     * This is called when writing a session, which can happen on the event thread, on the autosave timer thread, and
     * from the shutdown hook, so it loads the IGV genome list only -- never the much larger UCSC GenArk list -- and
     * reports a failure to the log rather than to a dialog.  If the list cannot be read the answer is false, and the
     * session records the expanded genome definition instead of the ID.
     *
     * @param genomeId
     * @return
     */
    public static synchronized boolean isIGVHosted(String genomeId) {
        loadIGVHostedIds();
        return genomeId != null && igvHostedIds.contains(genomeId);
    }

    /**
     * Populate the set of IGV hosted genome IDs from the IGV genome list, if that has not already happened as part of
     * reading the full record set.  A failure is logged and not retried.
     */
    private static void loadIGVHostedIds() {

        if (igvHostedIdsLoaded) {
            return;
        }
        igvHostedIdsLoaded = true;

        final IGVPreferences preferences = PreferencesManager.getPreferences();
        List<String> errors = new ArrayList<>();
        List<GenomeListItem> items = fetchGenomeList(preferences.get(GENOMES_SERVER_URL), "assembly", errors);
        if (items == null) {
            items = fetchGenomeList(preferences.get(BACKUP_GENOMES_SERVER_URL), "assembly", errors);
        }
        if (items != null) {
            recordIGVHostedIds(items);
        }
    }

    private static void recordIGVHostedIds(List<GenomeListItem> items) {
        for (GenomeListItem item : items) {
            if (item.getId() != null) {
                igvHostedIds.add(item.getId());
            }
        }
        igvHostedIdsLoaded = true;
    }

    public static synchronized GenomeListItem getGenomeListItem(String genomeId) {
        if (hostedGenomesMap == null) {
            hostedGenomesMap = new HashMap<>();
            for (GenomeListItem record : getRecords()) {
                hostedGenomesMap.put(record.getId(), record);
            }
        }
        return hostedGenomesMap.get(genomeId);
    }

private static List<GenomeListItem> readRecords() {

    records = new CopyOnWriteArrayList<>();

    final IGVPreferences preferences = PreferencesManager.getPreferences();
    final String genomesServerURL = preferences.get(GENOMES_SERVER_URL);
    final String backupGenomesServerURL = preferences.get(BACKUP_GENOMES_SERVER_URL);
    final String genarkURL = "https://hgdownload.soe.ucsc.edu/hubs/UCSC_GI.assemblyHubList.txt";

    List<String> errors = new ArrayList<>();

    // IGV hosted genome list
    List<GenomeListItem> igvGenomes = fetchGenomeList(genomesServerURL, "assembly", errors);
    if (igvGenomes == null) {
        log.error("Error loading genome list from: " + genomesServerURL);
        // Try backup server
        igvGenomes = fetchGenomeList(backupGenomesServerURL, "assembly", errors);
    }
    if (igvGenomes != null) {
        records.addAll(igvGenomes);
        recordIGVHostedIds(igvGenomes);
    }

    // UCSC Genark hosted genome list
    List<GenomeListItem> genarkGenomes = fetchGenomeList(genarkURL, "assembly", errors);
    if (genarkGenomes == null) {
        log.error("Error connecting to UCSC Genark server URL: " + genarkURL);
    } else {
        records.addAll(genarkGenomes);
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

/**
 * Fetch and parse a genome list, returning its records, or null if it could not be read.
 */
private static List<GenomeListItem> fetchGenomeList(String url, String idColumn, List<String> errors) {
    try {
        String genomeListContent = HttpUtils.getInstance().getContentsAsString(new URL(url));
        List<String> genomeListLines = Arrays.asList(genomeListContent.split("\\r?\\n"));
        String[] headers = parseHeaders(genomeListLines);
        return parseRecords(genomeListLines, headers, idColumn);
    } catch (Exception e) {
        log.error("Error loading genome list from: " + url, e);
        errors.add("Error loading genome list from: " + url + "   (" + e.getMessage() + ")");
        return null;
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


    private static List<GenomeListItem> parseRecords(List<String> genomeListLines, String [] headers, String idColumn) {

        List<GenomeListItem> items = new ArrayList<>();
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
                items.add(new GenomeListItem(displayableName, path, id, attributes));
            }
        }
        return items;
    }
}
