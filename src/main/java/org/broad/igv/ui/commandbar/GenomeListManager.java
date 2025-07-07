package org.broad.igv.ui.commandbar;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.event.GenomeResetEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.genome.GenomeDescriptor;
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

    private static final GenomeDescriptor DEFAULT_GENOME = new GenomeDescriptor(
            "Human (hg38)",
            "https://raw.githubusercontent.com/igvteam/igv-data/refs/heads/main/genomes/legacy/json/hg38.json",
            "hg38");

    private Map<String, GenomeDescriptor> genomeItemMap;

    private Map<String, GenomeDescriptor> remoteGenomesMap;

    private Map<String, GenomeDescriptor> downloadedGenomesMap;

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
    public Map<String, GenomeDescriptor> getGenomeItemMap() throws IOException {
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

    public List<GenomeDescriptor> getGenomeTableRecords() {
        List<GenomeDescriptor> items = new ArrayList<>(genomeItemMap.values());
        items.sort(sorter);
        return items;
    }


    /**
     * Add an item to the selectable genomes map and record in preferences.
     *
     * @param GenomeTableRecord
     */
    public void addGenomeItem(GenomeDescriptor GenomeTableRecord) {

        if (genomeItemMap.values().stream().anyMatch(item -> GenomeTableRecord.equals(item))) {
            return;
        }

        genomeItemMap.put(GenomeTableRecord.getId(), GenomeTableRecord);

        if (FileUtils.isRemote(GenomeTableRecord.getPath()) ||
                !DirectoryManager.getGenomeCacheDirectory().
                        equals(new File(GenomeTableRecord.getPath()).getParentFile())) {
            if (remoteGenomesMap == null) {
                remoteGenomesMap = new HashMap<>();
            }
            remoteGenomesMap.put(GenomeTableRecord.getId(), GenomeTableRecord);
            exportRemoteGenomesList();
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
    public GenomeDescriptor getGenomeTableRecord(String genomeId) {

        GenomeDescriptor matchingItem = genomeItemMap.get(genomeId);
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
     * If both .genome and .json files are found for the same genome ID the .json file is preferred, and the .genome
     * file deleted.   This complicates loading but is neccessary for a transition period.
     *
     * @return LinkedHashSet<GenomeTableRecord>
     * @throws IOException
     * @see GenomeDescriptor
     */
    public Map<String, GenomeDescriptor> getDownloadedGenomeMap() {

        if (downloadedGenomesMap == null) {

            downloadedGenomesMap = new HashMap<>();
            if (!DirectoryManager.getGenomeCacheDirectory().exists()) {
                return downloadedGenomesMap;
            }

            File[] files = DirectoryManager.getGenomeCacheDirectory().listFiles();

            // First loop - .genome files
            for (File file : files) {
                if (file.isDirectory()) {
                    continue;
                }
                if (file.getName().toLowerCase().endsWith(Globals.GENOME_FILE_EXTENSION)) {
                    ZipFile zipFile = null;
                    FileInputStream fis = null;
                    ZipInputStream zipInputStream = null;
                    try {
                        zipFile = new ZipFile(file);
                        fis = new FileInputStream(file);
                        zipInputStream = new ZipInputStream(new BufferedInputStream(fis));

                        ZipEntry zipEntry = zipFile.getEntry(org.broad.igv.feature.genome.load.GenomeDescriptor.GENOME_ARCHIVE_PROPERTY_FILE_NAME);
                        if (zipEntry == null) {
                            continue;    // Should never happen
                        }

                        InputStream inputStream = zipFile.getInputStream(zipEntry);
                        Properties properties = new Properties();
                        properties.load(inputStream);

                        GenomeDescriptor item =
                                new GenomeDescriptor(properties.getProperty(org.broad.igv.feature.genome.load.GenomeDescriptor.GENOME_ARCHIVE_NAME_KEY),
                                        file.getAbsolutePath(),
                                        properties.getProperty(org.broad.igv.feature.genome.load.GenomeDescriptor.GENOME_ARCHIVE_ID_KEY));

                        downloadedGenomesMap.put(item.getId(), item);

                    } catch (ZipException ex) {
                        log.error("\nZip error unzipping cached genome.", ex);
                        try {
                            file.delete();
                            zipInputStream.close();
                        } catch (Exception e) {
                            //ignore exception when trying to delete file
                        }
                    } catch (IOException ex) {
                        log.warn("\nIO error unzipping cached genome.", ex);
                        try {
                            file.delete();
                        } catch (Exception e) {
                            //ignore exception when trying to delete file
                        }
                    } finally {
                        try {
                            if (zipInputStream != null) {
                                zipInputStream.close();
                            }
                            if (zipFile != null) {
                                zipFile.close();
                            }
                            if (fis != null) {
                                fis.close();
                            }
                        } catch (IOException ex) {
                            log.warn("Error closing genome zip stream", ex);
                        }
                    }
                }
            }

            // Second loop - .json files
            for (File file : files) {
                if (file.isDirectory()) {
                    continue;
                }
                if (file.getName().toLowerCase().endsWith(".json")) {
                    try {
                        BufferedReader reader = new BufferedReader(new FileReader(file));
                        JsonParser parser = new JsonParser();
                        JsonElement rootElement = parser.parse(reader);
                        if (rootElement.isJsonObject()) {
                            JsonObject json = rootElement.getAsJsonObject();
                            JsonElement id = json.get("id");
                            JsonElement name = json.get("name");
                            JsonElement fasta = json.get("fastaURL");
                            JsonElement twobit = json.get("twoBitURL");
                            if (id == null) {
                                log.error("Error parsing " + file.getName() + ". \"id\" is required");
                                continue;
                            }
                            if (name == null) {
                                log.error("Error parsing " + file.getName() + ". \"name\" is required");
                                continue;
                            }
                            if (fasta == null && twobit == null) {
                                log.error("Error parsing " + file.getName() + ". One of either \"fastaURL\" or \"twoBitURL\" is required");
                                continue;
                            }

                            GenomeDescriptor item = new GenomeDescriptor(name.getAsString(), file.getAbsolutePath(), id.getAsString());
                            downloadedGenomesMap.put(item.getId(), item);

                        }
                    } catch (Exception e) {
                        log.error("Error parsing genome json: " + file.getAbsolutePath(), e);
                    }
                }
            }
        }

        return downloadedGenomesMap;
    }

    public void removeItems(List<GenomeDescriptor> removedValuesList) {

        boolean updateImportFile = false;
        for (GenomeDescriptor GenomeTableRecord : removedValuesList) {
            final String id = GenomeTableRecord.getId();
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
     * @see GenomeDescriptor
     */
    public Map<String, GenomeDescriptor> getRemoteGenomesMap() {

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
                            GenomeDescriptor item = new GenomeDescriptor(fields[0], file, fields[2]);
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
            for (GenomeDescriptor GenomeTableRecord : remoteGenomesMap.values()) {
                writer.print(GenomeTableRecord.getDisplayableName());
                writer.print("\t");
                writer.print(GenomeTableRecord.getPath());
                writer.print("\t");
                writer.println(GenomeTableRecord.getId());
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


    private static class GenomeListSorter implements Comparator<GenomeDescriptor> {

        @Override
        public int compare(GenomeDescriptor o1, GenomeDescriptor o2) {
            return o1.getDisplayableName().toLowerCase().compareTo(o2.getDisplayableName().toLowerCase());
        }
    }


    // Tacking on a timestamp & random number to avoid file collisions with parallel testing JVMs.  Not guaranteed unique
    // but highly unlikely to be repeated.
    public static final String TEST_REMOTE_GENOMES_FILE = "test-remote-genomes_" +
            System.currentTimeMillis() + "_" + Math.random() + ".txt";


}
