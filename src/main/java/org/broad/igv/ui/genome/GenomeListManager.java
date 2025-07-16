package org.broad.igv.ui.genome;

import com.google.gson.Gson;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.event.GenomeResetEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.genome.GenomeDownloadUtils;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.HostedGenomes;
import org.broad.igv.feature.genome.load.*;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;

/**
 * Singleton class for managing the list of genomes presented in the command bar combo box.
 */

public class GenomeListManager {

    private static Logger log = LogManager.getLogger(GenomeListManager.class);

    private static GenomeListManager theInstance;

    private static final String REMOTE_GENOMES_FILE = "remote-genomes.txt";

    private static final GenomeListItem DEFAULT_GENOME = new GenomeListItem(
            "Human (hg38)",
            "https://raw.githubusercontent.com/igvteam/igv-data/refs/heads/main/genomes/legacy/json/hg38.json",
            "hg38");

    private Map<String, GenomeListItem> genomeItemMap;

    private Map<String, GenomeListItem> remoteGenomesMap;

    private Map<String, GenomeListItem> downloadedGenomesMap;

    private static GenomeListSorter sorter = new GenomeListSorter();

    public synchronized static GenomeListManager getInstance() {
        if (theInstance == null) {
            theInstance = new GenomeListManager();
        }
        return theInstance;
    }

    private GenomeListManager() {
        genomeItemMap = new HashMap<>();
    }

    /**
     * Return the map of currently selectable genomes (id -> GenomeTableRecord).  If the map is null
     * initialize it from cached and user-defined genomes.
     *
     * @return
     * @throws IOException
     */
    public Map<String, GenomeListItem> getGenomeItemMap() throws IOException {
        if (genomeItemMap.isEmpty()) {
            rebuildGenomeItemMaps();
        }
        return genomeItemMap;
    }

    public void resetForTests() throws IOException {
        remoteGenomesMap = null;
        new File(TEST_REMOTE_GENOMES_FILE).delete();
        rebuildGenomeItemMaps();
    }

    public List<GenomeListItem> getGenomeTableRecords() {
        List<GenomeListItem> items = new ArrayList<>(genomeItemMap.values());
        items.sort(sorter);
        return items;
    }


    /**
     * Add an item to the selectable genomes map and record in preferences.
     *
     * @param genomeListItem
     */
    public void addGenomeItem(GenomeListItem genomeListItem) {

        if (genomeItemMap.values().stream().anyMatch(item -> genomeListItem.equals(item))) {
            return;
        }

        genomeItemMap.put(genomeListItem.getId(), genomeListItem);

        GenomeListItem hostedListItem = HostedGenomes.getGenomeListItem(genomeListItem.getId());
        boolean isHosted = hostedListItem != null && hostedListItem.equals(genomeListItem);
        boolean isCached = DirectoryManager.getGenomeCacheDirectory().
                equals(new File(genomeListItem.getPath()).getParentFile());

        // If this genome is otherwise unknown (not in the hosted set or local cache) record it in the remote genomes map
        if (!isHosted && !isCached) {
            if (remoteGenomesMap == null) {
                remoteGenomesMap = new HashMap<>();
            }

            if (hostedListItem == null || !hostedListItem.equals(genomeListItem)) {

                remoteGenomesMap.put(genomeListItem.getId(), genomeListItem);
                remoteGenomesMap.put(genomeListItem.getId(), genomeListItem);
                exportRemoteGenomesList();
            }
        }

        IGVEventBus.getInstance().post(new GenomeResetEvent());
    }


    /**
     * Searches through currently loaded GenomeTableRecords and returns
     * that with a matching ID. If not found, searches server and
     * user defined lists
     *
     * @param genomeId
     * @return
     */
    public GenomeListItem getGenomeTableRecord(String genomeId) {

        GenomeListItem matchingItem = genomeItemMap.get(genomeId);
        if (matchingItem == null) {

            // If not found rebuild the item map
            try {
                rebuildGenomeItemMaps();
            } catch (IOException e) {
                log.error("Error rebuilding genome item map", e);
            }
            matchingItem = genomeItemMap.get(genomeId);
        }
        return matchingItem;
    }

    /**
     * Rebuild the genome drop down info.
     */
    private void rebuildGenomeItemMaps() throws IOException {
        // Rebuild the selectable genomes map.  The order is imporant as downloaded genomes take precedence.

        remoteGenomesMap = null;
        downloadedGenomesMap = null;
        genomeItemMap.clear();
        genomeItemMap.putAll(getRemoteGenomesMap());
        genomeItemMap.putAll(getDownloadedGenomeMap());
        if (genomeItemMap.isEmpty()) {
            genomeItemMap.put(DEFAULT_GENOME.getId(), DEFAULT_GENOME);
        }
    }

    /**
     * Build and return a genome list item map from all local (downloaded) genome files (.genome and .json files) in
     * the igv/genomes directory.
     * <p>
     * If both .genome and .json files are found for the same genome ID the .json file is preferred.  This complicates
     * loading but is neccessary for a transition period.
     *
     * @return LinkedHashSet<GenomeTableRecord>
     * @throws IOException
     * @see GenomeListItem
     */
    public Map<String, GenomeListItem> getDownloadedGenomeMap() {

        if (downloadedGenomesMap == null) {

            downloadedGenomesMap = new HashMap<>();
            if (!DirectoryManager.getGenomeCacheDirectory().exists()) {
                return downloadedGenomesMap;
            }

            File[] files = DirectoryManager.getGenomeCacheDirectory().listFiles();

            // First loop - .json and .gbk files
            for (File file : files) {
                if (file.isDirectory()) {
                    continue;
                }
                if (file.getName().toLowerCase().endsWith(".json")) {

                    try (BufferedReader reader = new BufferedReader(new FileReader(file))) {

                        GenomeConfig config = (new Gson()).fromJson(reader, GenomeConfig.class);
                        String id = config.id;
                        if (id == null) {
                            log.error("GenomeConfig ID is null for file: " + file.getAbsolutePath());
                            continue;
                        }
                        String name = config.getName();
                        GenomeListItem item = new GenomeListItem(name, file.getAbsolutePath(), id);

                        if (GenomeUtils.isDeprecated(config)) {
                            GenomeListItem hostedItem = HostedGenomes.getGenomeListItem(id);
                            if (hostedItem != null) {
                                item = hostedItem;
                            }
                        }
                        downloadedGenomesMap.put(id, item);
                    } catch (Exception e) {
                        log.error("Error parsing genome json: " + file.getAbsolutePath(), e);
                    }
                } else if (file.getName().toLowerCase().endsWith(".gbk")) {
                    try {
                        String id = (new GenbankParser(file.getAbsolutePath())).getAccession();
                        String name = id;
                        GenomeListItem item = new GenomeListItem(name, file.getAbsolutePath(), id);
                        downloadedGenomesMap.put(id, item);
                    } catch (IOException e) {
                        log.error("Error parsing Genbank file: " + file.getAbsolutePath(), e);
                    }
                }
            }

            // Second loop - .genome files
            for (File file : files) {
                if (file.isDirectory()) {
                    continue;
                }
                if (file.getName().toLowerCase().endsWith(Globals.GENOME_FILE_EXTENSION)) {

                    String path = file.getAbsolutePath();
                    try (ZipFile zipFile = new ZipFile(path)) {
                        ZipEntry zipEntry = zipFile.getEntry(GenomeDescriptor.GENOME_ARCHIVE_PROPERTY_FILE_NAME);
                        if (zipEntry == null) {
                            throw new IOException("Missing genome archive property file.");
                        }

                        try (InputStream inputStream = zipFile.getInputStream(zipEntry)) {
                            Properties properties = new Properties();
                            properties.load(inputStream);

                            String id = properties.getProperty(GenomeDescriptor.GENOME_ARCHIVE_ID_KEY);
                            if (downloadedGenomesMap.containsKey(id)) {
                                log.info("Ignoring deprecated .genome file for genome: " + id);
                                continue;
                            }

                            String name = properties.getProperty(GenomeDescriptor.GENOME_ARCHIVE_NAME_KEY);
                            GenomeListItem genomeListItem = new GenomeListItem(name, path, id);

                            GenomeListItem hostedItem = HostedGenomes.getGenomeListItem(id);
                            String fastaURL = properties.getProperty(GenomeDescriptor.GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY);
                            if (hostedItem != null && fastaURL != null && GenomeUtils.isDeprecated(fastaURL)) {
                                log.warn("Genome file " + file.getName() + " is deprecated. Using hosted genome instead.");
                                genomeListItem = GenomeUtils.updateGenome(hostedItem);
                            }

                            downloadedGenomesMap.put(id, genomeListItem);
                        }
                    } catch (Exception e) {
                        log.warn("Error reading genome archive file: " + path + ". It may not be a valid genome archive.");
                    }
                }
            }
        }

        return downloadedGenomesMap;
    }


    public void removeItems(List<GenomeListItem> removedValuesList) {

        boolean updateImportFile = false;
        for (GenomeListItem genomeListItem : removedValuesList) {
            final String id = genomeListItem.getId();
            genomeItemMap.remove(id);
            if (remoteGenomesMap != null && remoteGenomesMap.containsKey(id)) {
                remoteGenomesMap.remove(id);
                updateImportFile = true;
            }
        }
        if (updateImportFile) {
            exportRemoteGenomesList();
        }
    }

    public void removeRemoteItem(String id) {
        if (remoteGenomesMap.containsKey(id)) {
            remoteGenomesMap.remove(id);
            exportRemoteGenomesList();
        }
    }

    private static BufferedReader getGenomeServerListReader() throws IOException {
        try {
            return ParsingUtils.openBufferedReader(
                    PreferencesManager.getPreferences().getGenomeListURL());
        } catch (IOException e) {
            log.error("Error fetching genome list: ", e);
            return ParsingUtils.openBufferedReader(
                    PreferencesManager.getPreferences().getBackupGenomeListURL());

        }
    }

    /**
     * Return the list of remote genomes for the dropdown.
     *
     * Remote genomes include all user-loaded genomes not stored locally in the igv/genomes folder.  This
     * does not include hosted genomes.
     *
     * @return LinkedHashSet<GenomeTableRecord>
     * @throws IOException
     * @see GenomeListItem
     */
    public Map<String, GenomeListItem> getRemoteGenomesMap() {

        if (remoteGenomesMap == null) {

            boolean updateClientGenomeListFile = false;

            remoteGenomesMap = new HashMap<>();

            File listFile = new File(DirectoryManager.getGenomeCacheDirectory(), getRemoteGenomesFilename());

            if (!listFile.exists()) {
                // Try old filename
                listFile = new File(DirectoryManager.getGenomeCacheDirectory(), "user-defined-genomes.txt");
            }

            if (listFile.exists()) {

                BufferedReader reader = null;

                boolean mightBeProperties = false;
                try {
                    reader = new BufferedReader(new FileReader(listFile));
                    String nextLine;
                    while ((nextLine = reader.readLine()) != null) {
                        if (nextLine.startsWith("#") || nextLine.trim().length() == 0) {
                            mightBeProperties = true;
                            continue;
                        }

                        // TODO: File should be called tsv maybe???
                        String[] fields = nextLine.split("\t");
                        if (fields.length < 3) {
                            if (mightBeProperties && fields[0].contains("=")) {
                                fields = nextLine.split("\\\\t");
                                if (fields.length < 3) {
                                    continue;
                                }
                                int idx = fields[0].indexOf("=");
                                fields[0] = fields[0].substring(idx + 1);
                            }
                        }

                        String file = fields[1];
                        // If a local file, see if it exists
                        if (!FileUtils.isRemote(file) && !FileUtils.resourceExists(file)) {
                            updateClientGenomeListFile = true;
                            continue;
                        }

                        try {
                            GenomeListItem item = new GenomeListItem(fields[0], file, fields[2]);
                            remoteGenomesMap.put(item.getId(), item);
                        } catch (Exception e) {
                            log.error("Error updating user genome list line '" + nextLine + "'", e);
                        }
                    }
                } catch (FileNotFoundException e) {
                    //We swallow this because the user may not have the file,
                    //which doesn't really matter
                    log.error(e);
                } catch (IOException e) {
                    log.error(e);
                    throw new RuntimeException(e);
                } finally {
                    if (reader != null) try {
                        reader.close();
                    } catch (IOException e) {
                    }
                }
                if (updateClientGenomeListFile) {
                    exportRemoteGenomesList();
                }
            }
        }
        return remoteGenomesMap;
    }


    /**
     * Export the list of selected remote genomes (genomes not in local cache)
     *
     * @throws IOException
     */
    public void exportRemoteGenomesList() {

        if (remoteGenomesMap == null) {
            return;
        }

        File listFile = new File(DirectoryManager.getGenomeCacheDirectory(), getRemoteGenomesFilename());
        File backup = null;
        if (listFile.exists()) {
            backup = new File(listFile.getAbsolutePath() + ".bak");
            try {
                FileUtils.copyFile(listFile, backup);
            } catch (IOException e) {
                log.error("Error backing up user-defined genome list file", e);
                backup = null;
            }
        }

        PrintWriter writer = null;
        try {
            writer = new PrintWriter(new BufferedWriter(new FileWriter(listFile)));
            for (GenomeListItem genomeListItem : remoteGenomesMap.values()) {
                writer.print(genomeListItem.getDisplayableName());
                writer.print("\t");
                writer.print(genomeListItem.getPath());
                writer.print("\t");
                writer.println(genomeListItem.getId());
            }

        } catch (Exception e) {
            if (backup != null) {
                try {
                    FileUtils.copyFile(backup, listFile);
                } catch (IOException e1) {
                    log.error("Error restoring genome-list file from backup");
                }
            }
            MessageUtils.showErrorMessage("Error updating user-defined genome list " + e.getMessage(), e);

        } finally {
            if (writer != null) writer.close();
            if (backup != null) backup.delete();
        }
    }


    private String getRemoteGenomesFilename() {
        if (Globals.isTesting()) {
            return TEST_REMOTE_GENOMES_FILE;
        } else {
            return REMOTE_GENOMES_FILE;
        }
    }


    private static class GenomeListSorter implements Comparator<GenomeListItem> {

        @Override
        public int compare(GenomeListItem o1, GenomeListItem o2) {
            return o1.getDisplayableName().toLowerCase().compareTo(o2.getDisplayableName().toLowerCase());
        }
    }


    // Tacking on a timestamp & random number to avoid file collisions with parallel testing JVMs.  Not guaranteed unique
// but highly unlikely to be repeated.
    public static final String TEST_REMOTE_GENOMES_FILE = "test-remote-genomes_" +
            System.currentTimeMillis() + "_" + Math.random() + ".txt";


}
