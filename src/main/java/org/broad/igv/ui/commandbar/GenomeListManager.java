package org.broad.igv.ui.commandbar;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGVMenuBar;
import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.Utilities;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;

/**
 * Created by jrobinso on 7/6/17.
 */
public class GenomeListManager {

    final static String GENOME_ARCHIVE_VERSION_KEY = "version";
    final static String GENOME_ARCHIVE_PROPERTY_FILE_NAME = "property.txt";
    final static String GENOME_ARCHIVE_ID_KEY = "id";
    final static String GENOME_ARCHIVE_NAME_KEY = "name";

    private static Logger log = LogManager.getLogger(GenomeManager.class);

    private static final String ACT_USER_DEFINED_GENOME_LIST_FILE = "user-defined-genomes.txt";
    public static final String TEST_USER_DEFINED_GENOME_LIST_FILE = GenomeManager.TEST_USER_DEFINED_GENOME_LIST_FILE;
    public static final GenomeListItem DEFAULT_GENOME = new GenomeListItem("Human (hg19)", "http://s3.amazonaws.com/igv.broadinstitute.org/genomes/hg19.genome", "hg19");

    private static GenomeListManager theInstance;

    private Map<String, GenomeListItem> genomeItemMap;

    private Map<String, GenomeListItem> userDefinedGenomeMap;

    private Map<String, GenomeListItem> serverGenomeMap;

    boolean serverGenomeListUnreachable = false;

    private static GenomeListSorter sorter = new GenomeListSorter();


    //private Map<String, GenomeListItem> genomeList;

    public synchronized static GenomeListManager getInstance() {
        if (theInstance == null) {
            theInstance = new GenomeListManager();
        }
        return theInstance;
    }

    private GenomeListManager() {
        genomeItemMap = new HashMap<>();
    }


////////////////////////////////////////////////////


    /**
     * Return the map of currently selectable genomes (id -> genomeListItem).  If the map is null
     * initialize it from cached and user-defined genomes.
     *
     * @return
     * @throws IOException
     */
    public Map<String, GenomeListItem> getGenomeItemMap() throws IOException {

        if (genomeItemMap.isEmpty()) {
            genomeItemMap.putAll(getUserDefinedGenomeMap());
            genomeItemMap.putAll(getCachedGenomeArchiveList());
            if (genomeItemMap.isEmpty()) {
                genomeItemMap.put(DEFAULT_GENOME.getId(), DEFAULT_GENOME);
            }
        }
        return genomeItemMap;
    }


    public List<GenomeListItem> getGenomeListItems() {
        List<GenomeListItem> items = new ArrayList<>(genomeItemMap.values());
        items.sort(sorter);
        return items;
    }


    /**
     * Completely rebuild the genome drop down info.
     */
    public void rebuildGenomeItemMap() throws IOException {
        genomeItemMap.clear();
        serverGenomeMap = null;
        userDefinedGenomeMap = null;
        getGenomeItemMap();
    }


    /**
     * Add an item to the selectable genomes map.  If not from server update the user defined file.
     *
     * @param genomeListItem
     * @param userDefined
     */
    public void addGenomeItem(GenomeListItem genomeListItem, boolean userDefined) {
        genomeItemMap.put(genomeListItem.getId(), genomeListItem);
        if (userDefined) {
            if (userDefinedGenomeMap == null) userDefinedGenomeMap = new HashMap<>();
            userDefinedGenomeMap.put(genomeListItem.getId(), genomeListItem);
            exportUserDefinedGenomeList();
        }

    }

    /**
     * Add a server-hosted genome list to the selectables map
     */
    public void addServerGenomeItem(GenomeListItem genomeListItem) {
        addGenomeItem(genomeListItem, false);
    }


    /**
     * Build a GenomeListItem from a path to a fasta file
     *
     * @param path
     * @return GenomeListItem representing this path, or null if file does not exist
     */
    public static GenomeListItem buildItemFromPath(String path) {

        String id = path;
        String name;
        if (HttpUtils.isRemoteURL(path)) {
            name = Utilities.getFileNameFromURL(path);
        } else {
            File file = new File(path);
            if (!file.exists()) {
                return null;
            }
            name = file.getName();
        }

        return new GenomeListItem(name, path, id);
    }

    /**
     * Return genome list item from currently selectable set. To search through
     * all server and user defined genomes, use {@link #getGenomeListItem(String)}
     *
     * @param genomeId
     * @return
     */
    public GenomeListItem getLoadedGenomeListItemById(String genomeId) {
        return genomeItemMap.get(genomeId);
    }

    /**
     * Searches through currently loaded GenomeListItems and returns
     * that with a matching ID. If not found, searches server and
     * user defined lists
     *
     * @param genomeId
     * @return
     */
    public GenomeListItem getGenomeListItem(String genomeId) {


        GenomeListItem matchingItem = genomeItemMap.get(genomeId);

        if (matchingItem == null || (System.currentTimeMillis() - matchingItem.getLastModified() >  GenomeManager.ONE_WEEK)) {

            // If genome archive was not found, or is more than 1 week old, check things not currently loaded
            matchingItem = getServerGenomeMap().get(genomeId);
            if (matchingItem != null) {
                return matchingItem;
            }

            // If still not found rebuild the item map
            try {
                rebuildGenomeItemMap();
            } catch (IOException e) {
                log.error("Error rebuilding genome item map", e);
            }
            matchingItem = genomeItemMap.get(genomeId);
        }
        return matchingItem;
    }

    /**
     * Gets a list of all the locally cached genome archive files.
     *
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    private Map<String, GenomeListItem> getCachedGenomeArchiveList()
            throws IOException {

        Map<String, GenomeListItem> cachedGenomeArchiveList = new HashMap<>();

        if (!DirectoryManager.getGenomeCacheDirectory().exists()) {
            return cachedGenomeArchiveList;
        }

        File[] files = DirectoryManager.getGenomeCacheDirectory().listFiles();
        for (File file : files) {

            if (file.isDirectory()) {
                continue;
            }

            if (!file.getName().toLowerCase().endsWith(Globals.GENOME_FILE_EXTENSION)) {
                continue;
            }


            ZipFile zipFile = null;
            FileInputStream fis = null;
            ZipInputStream zipInputStream = null;
            try {

                zipFile = new ZipFile(file);
                fis = new FileInputStream(file);
                zipInputStream = new ZipInputStream(new BufferedInputStream(fis));

                ZipEntry zipEntry = zipFile.getEntry(GENOME_ARCHIVE_PROPERTY_FILE_NAME);
                if (zipEntry == null) {
                    continue;    // Should never happen
                }

                InputStream inputStream = zipFile.getInputStream(zipEntry);
                Properties properties = new Properties();
                properties.load(inputStream);

                int version = 0;
                if (properties.containsKey(GENOME_ARCHIVE_VERSION_KEY)) {
                    try {
                        version = Integer.parseInt(
                                properties.getProperty(GENOME_ARCHIVE_VERSION_KEY));
                    } catch (Exception e) {
                        log.error("Error parsing genome version: " + version, e);
                    }
                }

                GenomeListItem item =
                        new GenomeListItem(properties.getProperty(GENOME_ARCHIVE_NAME_KEY),
                                file.getAbsolutePath(),
                                properties.getProperty(GENOME_ARCHIVE_ID_KEY));

                long lastModified = file.lastModified();
                item.setLastModified(lastModified);

                cachedGenomeArchiveList.put(item.getId(), item);
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
                    log.warn("Error closing genome zip stream!", ex);
                }
            }
        }

        return cachedGenomeArchiveList;
    }


    public Set<String> getServerGenomeIDs() {
        return getServerGenomeMap().keySet();
    }


    /**
     * Gets the collection of genome list items ids currently in use.
     *
     * @return Set of ids.
     */
    public Collection<String> getSelectableGenomeIDs() {
        return genomeItemMap.keySet();
    }


    public void removeAllItems(List<GenomeListItem> removedValuesList) {

        boolean updateImportFile = false;
        for (GenomeListItem genomeListItem : removedValuesList) {
            final String id = genomeListItem.getId();
            genomeItemMap.remove(id);
            if (userDefinedGenomeMap != null && userDefinedGenomeMap.containsKey(id)) {
                userDefinedGenomeMap.remove(id);
                updateImportFile = true;
            }
        }
        if (updateImportFile) {
            exportUserDefinedGenomeList();
        }
    }


    public void removeGenomeListItem(GenomeListItem genomeListItem) {

        final String id = genomeListItem.getId();
        genomeItemMap.remove(id);
        if (userDefinedGenomeMap != null && userDefinedGenomeMap.containsKey(id)) {
            userDefinedGenomeMap.remove(id);
            exportUserDefinedGenomeList();
        }
    }


    public List<GenomeListItem> getServerGenomeList() {
        List<GenomeListItem> items = new ArrayList<>(getServerGenomeMap().values());
        items.sort(sorter);
        return items;
    }

    /**
     * Gets a list of all the server genome archive files that
     * IGV knows about.
     *
     * @return List<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public Map<String, GenomeListItem> getServerGenomeMap() {

        if (serverGenomeListUnreachable) {
            return Collections.emptyMap();
        }

        if (serverGenomeMap == null) {
            serverGenomeMap = new HashMap<>();
            BufferedReader dataReader = null;
            InputStream inputStream = null;
            String genomeListURLString = "";
            try {
                genomeListURLString = PreferencesManager.getPreferences().getGenomeListURL();
                URL serverGenomeURL = HttpUtils.createURL(genomeListURLString);

                if (HttpUtils.isRemoteURL(genomeListURLString)) {
                    inputStream = HttpUtils.getInstance().openConnectionStream(serverGenomeURL);
                } else {
                    File file = new File(genomeListURLString.startsWith("file:") ? serverGenomeURL.getFile() : genomeListURLString);
                    inputStream = new FileInputStream(file);
                }


                dataReader = new BufferedReader(new InputStreamReader(inputStream));

                String genomeRecord;

                while ((genomeRecord = dataReader.readLine()) != null) {

                    if (genomeRecord.startsWith("<") || genomeRecord.startsWith("(#")) {
                        continue;
                    }

                    if (genomeRecord != null) {
                        genomeRecord = genomeRecord.trim();

                        String[] fields = genomeRecord.split("\t");
                        if ((fields != null) && (fields.length >= 3)) {


                            String name = fields[0];
                            String url = fields[1];
                            String id = fields[2];
                            GenomeListItem item = new GenomeListItem(name, url, id);
                            serverGenomeMap.put(item.getId(), item);

                        } else {
                            log.error("Found invalid server genome list record: " + genomeRecord);
                        }
                    }
                }
            } catch (Exception e) {
                serverGenomeListUnreachable = true;
                serverGenomeMap = Collections.emptyMap();
                log.error("Error fetching genome list: ", e);
                ConfirmDialog.optionallyShowInfoDialog("Warning: could not connect to the genome server (" +
                                genomeListURLString + ").    Only locally defined genomes will be available.",
                        Constants.SHOW_GENOME_SERVER_WARNING);

            } finally {
                if (dataReader != null) {
                    try {
                        dataReader.close();
                    } catch (IOException e) {
                        log.error(e);
                    }
                }
                if (inputStream != null) {
                    try {
                        inputStream.close();
                    } catch (IOException e) {
                        log.error(e);
                    }
                }
            }
        }

        if (IGVMenuBar.getInstance() != null) {
            IGVMenuBar.getInstance().notifyGenomeServerReachable(!serverGenomeListUnreachable);
        }
        return serverGenomeMap;
    }


    /**
     * Gets a list of all the user-defined genome archive files that
     * IGV knows about.
     *
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public Map<String, GenomeListItem> getUserDefinedGenomeMap() {

        if (userDefinedGenomeMap == null) {

            boolean updateClientGenomeListFile = false;

            userDefinedGenomeMap = new HashMap<>();

            File listFile = new File(DirectoryManager.getGenomeCacheDirectory(), getUserDefinedGenomeListFile());

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
                        if (!FileUtils.resourceExists(file)) {
                            updateClientGenomeListFile = true;
                            continue;
                        }

                        try {
                            GenomeListItem item = new GenomeListItem(fields[0], file, fields[2]);
                            userDefinedGenomeMap.put(item.getId(), item);
                        } catch (Exception e) {
                            log.error("Error updating user genome list line '" + nextLine + "'", e);
                        }
                    }
                } catch (FileNotFoundException e) {
                    //We swallow this because the user may not have the file,
                    //which doesn't really matter
                    log.info(e);
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
                    exportUserDefinedGenomeList();
                }
            }
        }
        return userDefinedGenomeMap;
    }


    /**
     * Export the user-define genome property file.
     *
     * @throws IOException
     */
    public void exportUserDefinedGenomeList() {

        if (userDefinedGenomeMap == null) {
            return;
        }

        File listFile = new File(DirectoryManager.getGenomeCacheDirectory(), getUserDefinedGenomeListFile());
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
            for (GenomeListItem genomeListItem : userDefinedGenomeMap.values()) {
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


    private String getUserDefinedGenomeListFile() {
        if (Globals.isTesting()) {
            return TEST_USER_DEFINED_GENOME_LIST_FILE;
        } else {
            return ACT_USER_DEFINED_GENOME_LIST_FILE;
        }

    }

    /**
     * Added for unit tests
     */
    public void clearUserDefinedGenomes() {
        userDefinedGenomeMap = null;
        new File(TEST_USER_DEFINED_GENOME_LIST_FILE).delete();
    }

    private static class GenomeListSorter implements Comparator<GenomeListItem> {

        @Override
        public int compare(GenomeListItem o1, GenomeListItem o2) {
            return o1.getDisplayableName().toLowerCase().compareTo(o2.getDisplayableName().toLowerCase());
        }
    }
}
