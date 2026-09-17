package org.igv.feature.genome;

import org.igv.DirectoryManager;
import org.igv.Globals;
import org.igv.feature.genome.load.HubGenomeLoader;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.genome.GenomeListItem;
import org.igv.ui.util.MessageUtils;
import org.igv.util.HttpUtils;

import java.io.File;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.function.Function;

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
     * Genomes listed by the IGV genome server (the GENOMES_SERVER_URL preference, or its backup), keyed by ID, or
     * null until first loaded.  These genomes can be restored from their ID alone, so sessions reference them by ID
     * rather than by an expanded genome definition.  There are a few dozen of them, as against the ~50,000 of the
     * UCSC GenArk list, so this is consulted before the full record set wherever an IGV hosted genome is the likely
     * answer.
     * <p>
     * Guarded by igvGenomeLock rather than the class monitor, which is held across the GenArk download.
     */
    private static Map<String, GenomeListItem> igvHostedGenomes;

    private static final Object igvGenomeLock = new Object();

    /**
     * Local copy of the last IGV genome list read, in the genome directory.  A live list is fetched every session;
     * this is only a fallback, so that an unreachable server does not cost a session its genome IDs.
     */
    private static final String IGV_GENOMES_CACHE_FILE = "hosted-genomes.tsv";

    /**
     * Errors from loading the IGV genome list, reported when the full record set is built.  The list may be loaded
     * silently, when writing a session, long before anything is in a position to show a dialog.
     */
    private static final List<String> igvGenomeErrors = new ArrayList<>();

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
    public static boolean isIGVHosted(String genomeId) {
        return igvHostedGenomes().containsKey(genomeId);
    }

    /**
     * Return the IGV hosted genome list keyed by ID, fetching it on first use.  A failure leaves the map empty and
     * is not retried, so one unreachable server does not cost every later call another timeout.
     */
    private static Map<String, GenomeListItem> igvHostedGenomes() {

        synchronized (igvGenomeLock) {
            if (igvHostedGenomes == null) {
                igvHostedGenomes = loadIGVHostedGenomes();
            }
            return igvHostedGenomes;
        }
    }

    /**
     * Read the IGV genome list, from the genome server if it can be reached and from the local copy if it cannot.
     * A list read from the server replaces the local copy.
     */
    private static Map<String, GenomeListItem> loadIGVHostedGenomes() {

        final IGVPreferences preferences = PreferencesManager.getPreferences();

        String content = fetchGenomeListContent(preferences.get(GENOMES_SERVER_URL), igvGenomeErrors);
        if (content == null) {
            content = fetchGenomeListContent(preferences.get(BACKUP_GENOMES_SERVER_URL), igvGenomeErrors);
        }

        List<GenomeListItem> items = parseGenomeList(content, "assembly", null);
        if (items == null || items.isEmpty()) {
            // Neither server could be read, or returned nothing usable.  Fall back to the last list read.
            items = parseGenomeList(readCachedGenomeList(), "assembly", null);
        } else {
            cacheGenomeList(content);
        }

        Map<String, GenomeListItem> genomes = new LinkedHashMap<>();
        if (items != null) {
            for (GenomeListItem item : items) {
                if (item.getId() != null) {
                    genomes.put(item.getId(), item);
                }
            }
        }
        return genomes;
    }

    /**
     * Return the contents of the local copy of the IGV genome list, or null if there is none to read.
     */
    private static String readCachedGenomeList() {
        File cacheFile = new File(DirectoryManager.getGenomeCacheDirectory(), IGV_GENOMES_CACHE_FILE);
        if (!cacheFile.exists()) {
            return null;
        }
        try {
            log.info("Using cached genome list: " + cacheFile.getAbsolutePath());
            return Files.readString(cacheFile.toPath(), StandardCharsets.UTF_8);
        } catch (Exception e) {
            log.error("Error reading cached genome list: " + cacheFile.getAbsolutePath(), e);
            return null;
        }
    }

    /**
     * Save the genome list just read, for use in a session that cannot reach the server.  A failure to write is of
     * no consequence beyond losing that fallback, so it is logged and otherwise ignored.
     */
    private static void cacheGenomeList(String content) {
        File cacheFile = new File(DirectoryManager.getGenomeCacheDirectory(), IGV_GENOMES_CACHE_FILE);
        try {
            Files.writeString(cacheFile.toPath(), content, StandardCharsets.UTF_8);
        } catch (Exception e) {
            log.error("Error caching genome list: " + cacheFile.getAbsolutePath(), e);
        }
    }

    /**
     * Return the record for a hosted genome ID, or null if there is none.
     * <p>
     * The IGV genome list is consulted first.  Only if the ID is not one of ours is the full record set built, which
     * means fetching the ~50,000 entry UCSC GenArk list -- so restoring a session that names an IGV hosted genome by
     * ID does not pay for it.  The two lists have no IDs in common, so looking at ours first changes no answer.
     *
     * @param genomeId
     * @return
     */
    public static synchronized GenomeListItem getGenomeListItem(String genomeId) {

        GenomeListItem item = igvHostedGenomes().get(genomeId);
        if (item != null) {
            return item;
        }

        if (hostedGenomesMap == null) {
            hostedGenomesMap = new HashMap<>();
            for (GenomeListItem record : getRecords()) {
                hostedGenomesMap.put(record.getId(), record);
            }
        }
        return hostedGenomesMap.get(genomeId);
    }

private static List<GenomeListItem> readRecords() {

    final String genarkURL = "https://hgdownload.soe.ucsc.edu/hubs/UCSC_GI.assemblyHubList.txt";

    List<GenomeListItem> allRecords = new ArrayList<>(igvHostedGenomes().values());
    List<String> errors = new ArrayList<>(igvGenomeErrors);

    // UCSC Genark hosted genome list.  These records are keyed by accession, which is the ID a Genark genome takes
    // when loaded -- the "genome" property of its hub is the accession.  The "assembly" column is a name such as
    // "Loxafr3.0", which is not what a session, batch command, or the last genome preference will name.  The list
    // has no "url" column either, so the path is derived from the accession.
    List<GenomeListItem> genarkGenomes =
            fetchGenomeList(genarkURL, "accession", errors, HubGenomeLoader::convertToHubURL);
    if (genarkGenomes != null) {
        allRecords.addAll(genarkGenomes);
    }

    if (!errors.isEmpty()) {
        StringBuilder sb = new StringBuilder();
        for (String error : errors) {
            sb.append(error).append("\n");
        }
        MessageUtils.showMessage(sb.toString());
    }

    return allRecords;
}

/**
 * Fetch and parse a genome list, returning its records, or null if it could not be read.
 *
 * @param pathFunction derives a record's path from its ID, for a list with no "url" column.  May be null.
 */
private static List<GenomeListItem> fetchGenomeList(String url, String idColumn, List<String> errors,
                                                    Function<String, String> pathFunction) {
    return parseGenomeList(fetchGenomeListContent(url, errors), idColumn, pathFunction);
}

/**
 * Fetch the contents of a genome list, returning null if the url could not be read.
 */
private static String fetchGenomeListContent(String url, List<String> errors) {
    try {
        return HttpUtils.getInstance().getContentsAsString(new URL(url));
    } catch (Exception e) {
        log.error("Error loading genome list from: " + url, e);
        errors.add("Error loading genome list from: " + url + "   (" + e.getMessage() + ")");
        return null;
    }
}

/**
 * Parse the contents of a genome list, returning its records, or null for null or unparsable content.
 *
 * @param pathFunction derives a record's path from its ID, for a list with no "url" column.  May be null.
 */
private static List<GenomeListItem> parseGenomeList(String content, String idColumn,
                                                    Function<String, String> pathFunction) {
    if (content == null) {
        return null;
    }
    try {
        List<String> genomeListLines = Arrays.asList(content.split("\\r?\\n"));
        String[] headers = parseHeaders(genomeListLines);
        return parseRecords(genomeListLines, headers, idColumn, pathFunction);
    } catch (Exception e) {
        log.error("Error parsing genome list", e);
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


    private static List<GenomeListItem> parseRecords(List<String> genomeListLines, String [] headers, String idColumn,
                                                     Function<String, String> pathFunction) {

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
                if (path == null && pathFunction != null && id != null) {
                    path = pathFunction.apply(id);
                }
                items.add(new GenomeListItem(displayableName, path, id, attributes));
            }
        }
        return items;
    }
}
