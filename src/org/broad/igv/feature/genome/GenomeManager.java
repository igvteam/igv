/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
/*
 * GenomeManager.java
 *
 * Created on November 9, 2007, 9:12 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.CytoBandFileParser;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.UIConstants;

import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.IGVHttpUtils;
import org.broad.igv.util.Utilities;

import javax.swing.*;
import java.io.*;
import java.net.*;
import java.util.*;
import java.util.zip.*;

/**
 * @author jrobinso
 */
public class GenomeManager {

    /**
     * Field description
     */
    public static final String SEQUENCE_FILE_EXTENSION = ".txt";
    /**
     * The refresh frequence in seconds.  Cached genomes will be refreshed
     * from the server if they are older than this value.
     */
    public static final long GENOME_REFRESH_FREQ = 7 * 24 * 3600;
    private static Logger log = Logger.getLogger(GenomeManager.class);
    private static GenomeManager theInstance = null;
    private Map<String, GenomeDescriptor> genomeDescriptorMap;
    private Map<String, Genome> genomes;

    final public static String USER_DEFINED_GENOME_LIST_FILE = "user-defined-genomes.txt";
    private static GenomeDescriptor DEFAULT;
    private String genomeId;
    public Genome genome;


    /**
     * Creates a new instance of GenomeManager
     */
    private GenomeManager() {
        genomeDescriptorMap = new HashMap();
        genomes = new Hashtable();
    }

    /**
     * Get the shared instance of the GenomeManager.
     *
     * @return GenomeManager
     */
    public static synchronized GenomeManager getInstance() {
        if (theInstance == null) {
            theInstance = new GenomeManager();
        }
        return theInstance;
    }

    /**
     * Return the genome identified by the id (e.g. mm8, hg17, etc).
     *
     * @param id Genome id.
     * @return
     */
    public Genome getGenome(String id) {

        Genome genome = genomes.get(id);
        if (genome == null) {
            GenomeDescriptor genomeDescriptor = genomeDescriptorMap.get(id);
            if (genomeDescriptor == null) {
                return null;
            } else {
                genome = loadGenome(genomeDescriptor);
                genomes.put(id, genome);

            }
        }
        return genome;
    }

    private Genome loadGenome(GenomeDescriptor genomeDescriptor) {
        InputStream is = null;
        try {

            InputStream inputStream = genomeDescriptor.getCytoBandStream();
            if (inputStream == null) {
                return null;
            }

            BufferedReader reader;
            if (genomeDescriptor.isCytoBandFileGZipFormat()) {
                is = new GZIPInputStream(inputStream);
                reader = new BufferedReader(new InputStreamReader(is));
            } else {
                is = new BufferedInputStream(inputStream);
                reader = new BufferedReader(new InputStreamReader(is));
            }

            Genome genome = new Genome(genomeDescriptor.getId());
            LinkedHashMap<String, Chromosome> chromMap = CytoBandFileParser.loadData(reader);
            genome.setChromosomeMap(chromMap, genomeDescriptor.isChromosomesAreOrdered());

            InputStream aliasStream = genomeDescriptor.getChrAliasStream();
            if (aliasStream != null) {
                try {
                    reader = new BufferedReader(new InputStreamReader(aliasStream));
                    genome.loadChrAliases(reader);
                } catch (Exception e) {
                    // We don't want to bomb if the alias load fails.  Just log it and proceed.
                    log.error("Error loading chromosome alias table");
                }
            }

            // Do this last so that user defined aliases have preference.  TODO This order dependence is rather fragile.
            genome.loadUserDefinedAliases();

            return genome;

        } catch (IOException ex) {
            log.warn("Error loading the genome!", ex);
            throw new RuntimeException("Error loading genome: " + genomeDescriptor.getName());
        } finally {
            try {
                if (is != null) {
                    is.close();
                }
            } catch (IOException ex) {
                log.warn("Error closing zip stream!", ex);
            }
        }
    }


    public GenomeListItem loadGenomeFromLocalFile(File zipFile) {

        GenomeListItem genomeListItem = null;

        try {
            boolean isUserDefined = false;
            File parentFolder = zipFile.getParentFile();
            if (parentFolder == null || !parentFolder.equals(Globals.getGenomeCacheDirectory())) {
                isUserDefined = true;
            }

            // user imported genomes are not versioned
            GenomeDescriptor genomeDescriptor = createGenomeDescriptor(zipFile, isUserDefined);
            genomeDescriptorMap.put(genomeDescriptor.getId(), genomeDescriptor);
            genomeListItem = new GenomeListItem(genomeDescriptor.getName(),
                    genomeDescriptor.getSequenceLocation(), genomeDescriptor.getId(),
                    genomeDescriptor.getVersion(),
                    isUserDefined);


        } catch (Exception e) {

            // Log any error and keep going
            log.warn("Error loading genome archive: " + zipFile, e);
        }

        return genomeListItem;
    }


    /**
     * Determine if the specified genome is already loaded and a descriptor
     * has been created for it.
     *
     * @param id Genome id.
     * @return true if genome is loaded.
     */
    public boolean isGenomeLoaded(String id) {
        return (genomeDescriptorMap.get(id) != null);
    }


    /**
     * Locates a genome in the set of known genome archive locations -
     * then loads it. This method is used ONLY for command line loading
     * of genomes.
     *
     * @param genome The genome to load.
     */
    public void findGenomeAndLoad(String genome) throws IOException {

        LinkedHashSet<GenomeListItem> genomes = getAllGenomeArchives(null);
        for (GenomeListItem item : genomes) {
            if (item.getId().equalsIgnoreCase(genome)) {
                String url = item.getLocation();
                if (!isGenomeLoaded(item.getId())) {
                    loadGenome(url, item.isUserDefined(), null);
                }
                break;
            }
        }
    }


    /**
     * Loads a system IGV genome archive into IGV's local genome cache and
     * then request it be loaded. If the genome is user-define it is loaded
     * directly from it's own location.
     *
     * @param item          A GenomeListItem.
     * @param isUserDefined true if a user genome.
     * @throws FileNotFoundException
     * @see GenomeListItem
     */
    public void loadGenome(GenomeListItem item, boolean isUserDefined)
            throws FileNotFoundException {


        if (log.isDebugEnabled()) {
            log.debug("Enter loadGenome");
        }

        if (isGenomeLoaded(item.getId())) {
            return;
        }
        loadGenome(item.getLocation(), isUserDefined, null);
    }


    /**
     * Copies a system IGV genome archive into IGV's local genome cache and
     * then request it be loaded. If the genome is user-define it is loaded
     * directly from it's own location.
     *
     * @param genomeArchiveFileLocation
     * @param isUserDefined             true if a user genome.
     * @param monitor                   ProgressMonitor
     * @return GenomeListItem.
     * @throws FileNotFoundException
     * @see GenomeListItem
     */
    public GenomeListItem loadGenome(
            String genomeArchiveFileLocation,
            boolean isUserDefined,
            ProgressMonitor monitor)
            throws FileNotFoundException {

        if (log.isDebugEnabled()) {
            log.debug("Enter loadGenome");
        }


        GenomeListItem genomeListItem = null;

        if (genomeArchiveFileLocation == null) {
            return getDefaultGenomeListItem();
        }

        if (!genomeArchiveFileLocation.trim().endsWith(Globals.GENOME_FILE_EXTENSION)) {
            throw new RuntimeException(
                    "The extension of archive [" + genomeArchiveFileLocation + "] is not an IGV genome archive extension");
        }

        try {

            File archiveFile;

            // If a System genome copy to local cache
            if (!isUserDefined) {

                if (!Globals.getGenomeCacheDirectory().exists()) {
                    Globals.getGenomeCacheDirectory().mkdir();
                }

                if (IGVHttpUtils.isURL(genomeArchiveFileLocation.toLowerCase())) {

                    URL genomeArchiveURL = new URL(genomeArchiveFileLocation);
                    String fileName = Utilities.getFileNameFromURL(
                            URLDecoder.decode(new URL(genomeArchiveFileLocation).getFile(), "UTF-8"));

                    archiveFile = new File(Globals.getGenomeCacheDirectory(), fileName);

                    refreshCache(archiveFile, genomeArchiveURL);

                } else {

                    archiveFile = new File(genomeArchiveFileLocation);
                }

                if (monitor != null) {
                    monitor.fireProgressChange(25);
                }

                if (!archiveFile.exists()) {
                    throw new FileNotFoundException(genomeArchiveFileLocation);
                }


            } else {

                // New genome archive file
                archiveFile = new File(genomeArchiveFileLocation);

                if (!archiveFile.exists()) {
                    throw new FileNotFoundException(genomeArchiveFileLocation);
                }

                File userDefinedListFile = new File(Globals.getGenomeCacheDirectory(), USER_DEFINED_GENOME_LIST_FILE);
                Properties listProperties =
                        retrieveUserDefinedGenomeListFromFile(userDefinedListFile);

                if (listProperties == null) {
                    listProperties = new Properties();
                }

                String record = buildClientSideGenomeListRecord(archiveFile, isUserDefined);

                if (record != null) {
                    int version = 0;
                    String[] fields = record.split("\t");
                    listProperties.setProperty(fields[2], record);
                    GenomeImporter.storeUserDefinedGenomeListToFile(userDefinedListFile, listProperties);
                    genomeListItem = new GenomeListItem(fields[0], fields[1], fields[2], version, isUserDefined);
                }

            }

            // If archive file exists load  it into IGV

            if (monitor != null) {
                monitor.fireProgressChange(25);
            }

            loadGenomeFromLocalFile(archiveFile);

            if (monitor != null) {
                monitor.fireProgressChange(25);
            }

        } catch (MalformedURLException e) {
            log.warn("Error Importing Genome!", e);
        } catch (FileNotFoundException e) {
            throw e;
        } catch (SocketException e) {
            throw new GenomeServerException("Server connection error", e);
        } catch (IOException e) {
            log.warn("Error Importing Genome!", e);
        }

        if (log.isDebugEnabled()) {
            log.debug("Exit loadGenome");
        }

        return genomeListItem;
    }


    /**
     * Refresh a locally cached genome
     *
     * @param archiveFile
     * @param genomeArchiveURL
     * @throws IOException
     */
    private void refreshCache(File archiveFile, URL genomeArchiveURL) {
        // Look in cache first

        InputStream is = null;

        try {
            if (archiveFile.exists()) {
                // Force an restorePersistentState of cached genomes periodically
                boolean forceUpdate = ((System.currentTimeMillis() - archiveFile.lastModified()) / 1000) > GENOME_REFRESH_FREQ;
                if (forceUpdate) {
                    log.info("Refreshing genome: " + genomeArchiveURL.toString());
                    File tmpFile = new File(archiveFile.getAbsolutePath() + ".tmp");
                    is = IGVHttpUtils.openConnectionStream(genomeArchiveURL);
                    FileUtils.createFileFromStream(is, tmpFile);
                    FileUtils.copyFile(tmpFile, archiveFile);
                    tmpFile.deleteOnExit();
                }
            } else {
                // Copy file directly from the server to local cache.
                is = IGVHttpUtils.openConnectionStream(genomeArchiveURL);
                FileUtils.createFileFromStream(is, archiveFile);
            }
        }
        catch (Exception e) {
            log.error("Error refreshing genome cache. ", e);
            MessageUtils.showMessage(("An error was encountered refreshing the genome cache: " + e.getMessage() +
                    "<br> If this problem persists please contact igv-help@broadinstitute.org"));
        }
        finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    log.error("Error closing genome stream: " + genomeArchiveURL.toString(), e);
                }
            }
        }
    }


    /**
     * Creates a genome descriptor.
     */
    private GenomeDescriptor createGenomeDescriptor(File f, boolean userDefined)
            throws IOException {


        String zipFilePath = f.getAbsolutePath();

        if (!f.exists()) {
            log.error("Genome file: " + f.getAbsolutePath() + " does not exist.");
            return null;
        }


        GenomeDescriptor genomeDescriptor = null;
        Map<String, ZipEntry> zipEntries = new HashMap();
        ZipFile zipFile = new ZipFile(zipFilePath);


        ZipInputStream zipInputStream = null;
        try {
            zipInputStream = new ZipInputStream(new FileInputStream(f));
            ZipEntry zipEntry = zipInputStream.getNextEntry();

            while (zipEntry != null) {
                String zipEntryName = zipEntry.getName();
                zipEntries.put(zipEntryName, zipEntry);

                if (zipEntryName.equalsIgnoreCase(Globals.GENOME_ARCHIVE_PROPERTY_FILE_NAME)) {
                    InputStream inputStream = zipFile.getInputStream(zipEntry);
                    Properties properties = new Properties();
                    properties.load(inputStream);

                    // Cytoband
                    String cytobandZipEntryName = properties.getProperty(Globals.GENOME_ARCHIVE_CYTOBAND_FILE_KEY);

                    // RefFlat
                    String geneFileName = properties.getProperty(Globals.GENOME_ARCHIVE_GENE_FILE_KEY);

                    String chrAliasFileName = properties.getProperty(Globals.GENOME_CHR_ALIAS_FILE_KEY);

                    String sequenceLocation = properties.getProperty(Globals.GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY);

                    if ((sequenceLocation != null) && !IGVHttpUtils.isURL(sequenceLocation)) {
                        File tempZipFile = new File(zipFilePath);
                        File sequenceFolder = new File(tempZipFile.getParent(), sequenceLocation);
                        sequenceLocation = sequenceFolder.getCanonicalPath();
                        sequenceLocation.replace('\\', '/');
                    }


                    int version = 0;
                    String versionString = properties.getProperty(Globals.GENOME_ARCHIVE_VERSION_KEY);
                    if (versionString != null) {
                        try {
                            version = Integer.parseInt(versionString);
                        } catch (Exception e) {
                            log.error("Error parsing version string: " + versionString);
                        }
                    }

                    boolean chrNamesAltered = false;
                    String chrNamesAlteredString = properties.getProperty("filenamesAltered");
                    if (chrNamesAlteredString != null) {
                        try {
                            chrNamesAltered = Boolean.parseBoolean(chrNamesAlteredString);
                        } catch (Exception e) {
                            log.error("Error parsing version string: " + versionString);
                        }
                    }

                    boolean chromosomesAreOrdered = false;
                    String tmp = properties.getProperty(Globals.GENOME_ORDERED_KEY);
                    if (tmp != null) {
                        try {
                            chromosomesAreOrdered = Boolean.parseBoolean(tmp);
                        }
                        catch (Exception e) {
                            log.error("Error parsing ordered string: " + tmp);
                        }
                    }

                    // The new descriptor
                    genomeDescriptor = new GenomeZipDescriptor(
                            properties.getProperty(Globals.GENOME_ARCHIVE_NAME_KEY),
                            version,
                            chrNamesAltered,
                            properties.getProperty(Globals.GENOME_ARCHIVE_ID_KEY),
                            cytobandZipEntryName,
                            geneFileName,
                            chrAliasFileName,
                            properties.getProperty(Globals.GENOME_GENETRACK_NAME, "Gene"),
                            sequenceLocation,
                            zipFile,
                            zipEntries,
                            userDefined,
                            chromosomesAreOrdered);

                }
                zipEntry = zipInputStream.getNextEntry();
            }
        } finally {
            try {
                if (zipInputStream != null) {
                    zipInputStream.close();
                }
            } catch (IOException ex) {
                log.warn("Error closing imported genome zip stream!", ex);
            }
        }
        return genomeDescriptor;
    }


    /**
     * Gets the descriptor for a specific genome.
     *
     * @param id
     * @return GenomeDescriptor
     */
    public GenomeDescriptor getGenomeDescriptor(String id) {
        if (genomeDescriptorMap.containsKey(id)) {
            return genomeDescriptorMap.get(id);
        } else {
            return getDefaultGenomeDescriptor();
        }
    }


    boolean serverGenomeListUnreachable = false;

    /**
     * Gets a list of all the server genome archive files that
     * IGV knows about.
     *
     * @param excludedArchivesUrls The set of file location to exclude in the return list.
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public LinkedHashSet<GenomeListItem> getServerGenomeArchiveList(Set excludedArchivesUrls)
            throws IOException {

        if (serverGenomeListUnreachable) {
            return null;
        }

        LinkedHashSet<GenomeListItem> genomeItemList = new LinkedHashSet();
        BufferedReader dataReader = null;
        InputStream inputStream = null;
        String genomeListURL = "";
        try {
            genomeListURL = PreferenceManager.getInstance().getGenomeListURL();
            URL serverGenomeArchiveList = new URL(genomeListURL);

            if (IGVHttpUtils.isURL(genomeListURL)) {
                inputStream = IGVHttpUtils.openConnectionStream(serverGenomeArchiveList);
            } else {
                File file = new File(genomeListURL.startsWith("file:") ? serverGenomeArchiveList.getFile() : genomeListURL);
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

                        // Throw away records we don't want to see
                        if (excludedArchivesUrls != null) {
                            if (excludedArchivesUrls.contains(fields[1])) {
                                continue;
                            }
                        }
                        int version = 0;
                        if (fields.length > 3) {
                            try {
                                version = Integer.parseInt(fields[3]);
                            } catch (Exception e) {
                                log.error("Error parsing genome version: " + fields[0], e);
                            }
                        }

                        try {
                            GenomeListItem item = new GenomeListItem(fields[0], fields[1], fields[2], version, false);
                            genomeItemList.add(item);
                        } catch (Exception e) {
                            log.error(
                                    "Error reading a line from server genome list" + " line was: [" + genomeRecord + "]",
                                    e);
                        }
                    } else {
                        log.error("Found invalid server genome list record: " + genomeRecord);
                    }
                }
            }
        } catch (Exception e) {
            serverGenomeListUnreachable = true;
            log.error("Error fetching genome list: ", e);
            ConfirmDialog.optionallyShowInfoDialog("Warning: could not connect to the genome server (" +
                    genomeListURL + ").    Only locally defined genomes will be available.",
                    PreferenceManager.SHOW_GENOME_SERVER_WARNING);
        } finally {
            if (dataReader != null) {
                dataReader.close();
            }
            if (inputStream != null) {
                inputStream.close();
            }
        }
        return genomeItemList;
    }

    /**
     * Gets a list of all the user-defined genome archive files that
     * IGV knows about.
     *
     * @param excludedArchivesUrls The set of file location to exclude in the
     *                             return list.
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public LinkedHashSet<GenomeListItem> getUserDefinedGenomeArchiveList(Set excludedArchivesUrls)
            throws IOException {

        boolean clientGenomeListNeedsRebuilding = false;
        LinkedHashSet<GenomeListItem> genomeItemList = new LinkedHashSet();

        File listFile = new File(Globals.getGenomeCacheDirectory(), USER_DEFINED_GENOME_LIST_FILE);

        Properties listProperties = retrieveUserDefinedGenomeListFromFile(listFile);

        if (listProperties != null) {

            Collection records = listProperties.values();

            for (Object value : records) {

                String record = (String) value;
                if (record.trim().equals("")) {
                    continue;
                }

                String[] fields = record.split("\t");

                // Throw away records we don't want to see
                if (excludedArchivesUrls != null) {
                    if (excludedArchivesUrls.contains(fields[1])) {
                        continue;
                    }
                }

                File file = new File(fields[1]);
                if (file.isDirectory()) {
                    continue;
                }
                if (!file.exists()) {
                    clientGenomeListNeedsRebuilding = true;
                    continue;
                }

                if (!file.getName().toLowerCase().endsWith(Globals.GENOME_FILE_EXTENSION)) {
                    continue;
                }
                GenomeListItem item = new GenomeListItem(fields[0], file.getAbsolutePath(),
                        fields[2], 0, true);
                genomeItemList.add(item);
            }
        }
        if (clientGenomeListNeedsRebuilding) {
            rebuildClientGenomeList(genomeItemList);
        }
        return genomeItemList;
    }

    /**
     * Method description
     */
    public void clearGenomeCache() {

        File[] files = Globals.getGenomeCacheDirectory().listFiles();
        for (File file : files) {
            if (file.getName().toLowerCase().endsWith(Globals.GENOME_FILE_EXTENSION)) {
                file.delete();
            }
        }

    }

    /**
     * Gets a list of all the locally cached genome archive files that
     * IGV knows about.
     *
     * @param excludedArchivesUrls The set of file location to exclude in the
     *                             return list.
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public LinkedHashSet<GenomeListItem> getCachedGenomeArchiveList(Set excludedArchivesUrls)
            throws IOException {

        LinkedHashSet<GenomeListItem> genomeItemList = new LinkedHashSet();

        if (!Globals.getGenomeCacheDirectory().exists()) {
            return genomeItemList;
        }

        File[] files = Globals.getGenomeCacheDirectory().listFiles();
        for (File file : files) {

            if (file.isDirectory()) {
                continue;
            }

            if (!file.getName().toLowerCase().endsWith(Globals.GENOME_FILE_EXTENSION)) {
                continue;
            }

            URL zipUrl = file.toURI().toURL();

            // Throw away records we don't want to see
            if (excludedArchivesUrls != null) {
                if (excludedArchivesUrls.contains(URLDecoder.decode(zipUrl.getFile(), "UTF-8"))) {
                    continue;
                }
            }

            ZipFile zipFile = null;
            FileInputStream fis = null;
            ZipInputStream zipInputStream = null;
            try {

                zipFile = new ZipFile(file);
                fis = new FileInputStream(file);
                zipInputStream = new ZipInputStream(new BufferedInputStream(fis));

                ZipEntry zipEntry = zipFile.getEntry(Globals.GENOME_ARCHIVE_PROPERTY_FILE_NAME);
                if (zipEntry == null) {
                    continue;    // Should never happen
                }

                InputStream inputStream = zipFile.getInputStream(zipEntry);
                Properties properties = new Properties();
                properties.load(inputStream);

                int version = 0;
                if (properties.containsKey(Globals.GENOME_ARCHIVE_VERSION_KEY)) {
                    try {
                        version = Integer.parseInt(
                                properties.getProperty(Globals.GENOME_ARCHIVE_VERSION_KEY));
                    } catch (Exception e) {
                        log.error("Error parsing genome version: " + version, e);
                    }
                }

                GenomeListItem item =
                        new GenomeListItem(properties.getProperty(Globals.GENOME_ARCHIVE_NAME_KEY),
                                file.getAbsolutePath(),
                                properties.getProperty(Globals.GENOME_ARCHIVE_ID_KEY),
                                version,
                                false);
                genomeItemList.add(item);
            } catch (ZipException ex) {
                log.error("\nZip error unzipping cached genome.", ex);
                JOptionPane.showMessageDialog(null, "Fatal error loading genome file: " + file.getAbsolutePath() +
                        "\n     *** " + ex.getMessage() + " ***" +
                        "\nIGV must exit.  If this problem persists contact igv-help@broadinstitute.org");
                try {
                    file.deleteOnExit();
                    zipInputStream.close();
                }
                catch (Exception e) {
                    //ignore exception when trying to delete file
                }
                System.exit(-1);


            } catch (IOException ex) {
                log.warn("\nIO error unzipping cached genome.", ex);
                try {
                    file.delete();
                }
                catch (Exception e) {
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

        return genomeItemList;
    }

    /**
     * Gets a list of all the server and client-side genome archive files that
     * IGV knows about.
     *
     * @param excludedArchivesUrls The set of file location to exclude in the
     *                             return list.
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public LinkedHashSet<GenomeListItem> getAllGenomeArchives(Set excludedArchivesUrls)
            throws IOException {

        LinkedHashSet<GenomeListItem> genomeListItems = new LinkedHashSet();

        // Build a single available genome list from both client, server
        // and cached information. This allows us to process
        // everything the same way.
        LinkedHashSet<GenomeListItem> serverSideItemList = null;
        try {
            serverSideItemList = getServerGenomeArchiveList(null);
        } catch (UnknownHostException e) {
            log.error(UIConstants.CANNOT_ACCESS_SERVER_GENOME_LIST, e);
        } catch (SocketException e) {
            log.error(UIConstants.CANNOT_ACCESS_SERVER_GENOME_LIST, e);
        }

        LinkedHashSet<GenomeListItem> cacheGenomeItemList = getCachedGenomeArchiveList(excludedArchivesUrls);

        LinkedHashSet<GenomeListItem> userDefinedItemList = getUserDefinedGenomeArchiveList(excludedArchivesUrls);


        if (serverSideItemList != null) {
            genomeListItems.addAll(serverSideItemList);
        }
        if (userDefinedItemList != null) {
            genomeListItems.addAll(userDefinedItemList);
        }
        if (cacheGenomeItemList != null) {
            genomeListItems.addAll(cacheGenomeItemList);
        }

        if (genomeListItems.isEmpty()) {
            GenomeDescriptor defaultDes = getDefaultGenomeDescriptor();
            GenomeListItem defaultItem = getDefaultGenomeListItem();
            genomeListItems.add(defaultItem);
            Genome genome = loadGenome(defaultDes);
            if (genome != null) {
                genomes.put(defaultDes.getId(), genome);
            }
        }

        return genomeListItems;
    }

    /**
     * Reconstructs the user-define genome property file.
     *
     * @param genomeItemList The list of user-define genome GenomeListItem
     *                       objects to store in the property file.
     * @throws IOException
     */
    public void rebuildClientGenomeList(LinkedHashSet<GenomeListItem> genomeItemList)
            throws IOException {

        if ((genomeItemList == null)) {
            return;
        }

        File listFile = new File(Globals.getGenomeCacheDirectory(), USER_DEFINED_GENOME_LIST_FILE);

        if (!listFile.exists()) {
            listFile.createNewFile();
        }

        StringBuffer buffer = new StringBuffer();
        Properties listProperties = new Properties();
        for (GenomeListItem genomeListItem : genomeItemList) {

            buffer.append(genomeListItem.getDisplayableName());
            buffer.append("\t");
            buffer.append(genomeListItem.getLocation());
            buffer.append("\t");
            buffer.append(genomeListItem.getId());

            listProperties.setProperty(genomeListItem.getId(), buffer.toString());
            buffer.delete(0, buffer.length());
        }
        GenomeImporter.storeUserDefinedGenomeListToFile(listFile, listProperties);
    }

    /**
     * Create a genonem list record (same format as used by the server) for
     * genome files.
     *
     * @param genomeArchive The genome file from which to extract a record.
     * @param userDefined   true if archive is a user-defined genome archive.
     * @return A tab delimetered genome list record containing
     *         ( name[tab]genome location[tab]genomeId ).
     * @throws IOException
     */
    public String buildClientSideGenomeListRecord(File genomeArchive, boolean userDefined)
            throws IOException {

        GenomeDescriptor genomeDescriptor = createGenomeDescriptor(genomeArchive, userDefined);

        StringBuffer buffer = new StringBuffer();
        buffer.append(genomeDescriptor.getName());
        buffer.append("\t");
        buffer.append(genomeArchive.getAbsoluteFile());
        buffer.append("\t");
        buffer.append(genomeDescriptor.getId());
        return buffer.toString();
    }

    /**
     * Read the user-defined genome property file to find enough information to
     * display the genome in IGV.
     *
     * @param file A java properties file containing tab delimetered data
     *             (display name [tab] genome file location [tab] genome id) about
     *             the user-defined genome.
     * @return A java Properties object contain the file's content.
     */
    public static Properties retrieveUserDefinedGenomeListFromFile(File file) {

        Properties properties = new Properties();

        if ((file != null) && file.exists()) {
            FileInputStream input = null;
            try {
                input = new FileInputStream(file);
                properties.load(input);
            } catch (FileNotFoundException e) {
                log.error("Property file for user-defined genomes was not " + "found!", e);
            } catch (IOException e) {
                log.error("Error readin property file for user-defined " + "genomes!", e);
            } finally {
                if (input != null) {
                    try {
                        input.close();
                    } catch (IOException e) {
                        log.error("Error closing property file for " + "user-defined genomes!",
                                e);
                    }
                }
            }
        }
        return properties;
    }

    /**
     * Create an IGV representation of a user-defined genome.
     *
     * @param genomeZipLocation              A File path to a directory in which the .genome
     *                                       output file will be written.
     * @param cytobandFileName               A File path to a file that contains cytoband data.
     * @param refFlatFileName                A File path to a gene file.
     * @param fastaFileName                  A File path to a FASTA file, a .gz file containing a
     *                                       single FASTA file, or a directory containing ONLY FASTA files.
     * @param relativeSequenceLocation       A relative path to the location
     *                                       (relative to the .genome file to be created) where the sequence data for
     *                                       the new genome will be written.
     * @param genomeDisplayName              The unique user-readable name of the new genome.
     * @param genomeId                       The id to be assigned to the genome.
     * @param genomeFileName                 The file name (not path) of the .genome archive
     *                                       file to be created.
     * @param monitor                        A ProgressMonitor used to track progress - null,
     *                                       if no progress updating is required.
     * @param sequenceOutputLocationOverride
     * @return GenomeListItem
     * @throws FileNotFoundException
     */
    public GenomeListItem defineGenome(String genomeZipLocation,
                                       String cytobandFileName,
                                       String refFlatFileName,
                                       String fastaFileName,
                                       String chrAliasFileName,
                                       String relativeSequenceLocation,
                                       String genomeDisplayName,
                                       String genomeId,
                                       String genomeFileName,
                                       ProgressMonitor monitor,
                                       String sequenceOutputLocationOverride)
            throws IOException {

        File archiveFile;
        File zipFileLocation = null;
        File fastaInputFile = null;
        File refFlatFile = null;
        File cytobandFile = null;
        File chrAliasFile = null;
        File sequenceLocation;

        if ((genomeZipLocation != null) && (genomeZipLocation.trim().length() != 0)) {
            zipFileLocation = new File(genomeZipLocation);
            PreferenceManager.getInstance().setLastGenomeImportDirectory(zipFileLocation);
        }

        if ((cytobandFileName != null) && (cytobandFileName.trim().length() != 0)) {
            cytobandFile = new File(cytobandFileName);
        }

        if ((refFlatFileName != null) && (refFlatFileName.trim().length() != 0)) {
            refFlatFile = new File(refFlatFileName);
        }

        if ((chrAliasFileName != null) && (chrAliasFileName.trim().length() != 0)) {
            chrAliasFile = new File(chrAliasFileName);
        }

        if ((fastaFileName != null) && (fastaFileName.trim().length() != 0)) {
            fastaInputFile = new File(fastaFileName);

            // The sequence info only matters if we have FASTA
            if ((relativeSequenceLocation != null) && (relativeSequenceLocation.trim().length() != 0)) {
                sequenceLocation = new File(genomeZipLocation, relativeSequenceLocation);
                if (!sequenceLocation.exists()) {
                    sequenceLocation.mkdir();
                }
            }
        }

        if (monitor != null) {
            monitor.fireProgressChange(25);
        }

        archiveFile = (new GenomeImporter()).createGenomeArchive(zipFileLocation,
                genomeFileName, genomeId, genomeDisplayName, relativeSequenceLocation,
                fastaInputFile, refFlatFile, cytobandFile, chrAliasFile,
                sequenceOutputLocationOverride, monitor);

        if (monitor != null) {
            monitor.fireProgressChange(75);
        }


        if (log.isDebugEnabled()) {
            log.debug("Call loadGenome");
        }

        return loadGenome(archiveFile.getAbsolutePath(), true, null);
    }


    /**
     * This method takes any string and generates a unique key based on
     * that string. Other class typically call this method to make a genome id
     * from some string but it can be used to generate a key for anything.
     * TODO This method probably belongs in a Utilities class.
     *
     * @param text
     * @return A unique key based on the text.
     */
    public static String generateGenomeKeyFromText(String text) {

        if (text == null) {
            return null;
        }

        text = text.trim();
        text = text.replace(" ", "_");
        text = text.replace("(", "");
        text = text.replace(")", "");
        text = text.replace("[", "");
        text = text.replace("]", "");
        text = text.replace("{", "");
        text = text.replace("}", "");
        return text;
    }

    public GenomeDescriptor getDefaultGenomeDescriptor() {
        if (DEFAULT == null) {
            DEFAULT = new GenomeResourceDescriptor("Human hg18", 0, "hg18",
                    "/resources/hg18_cytoBand.txt", null, null, null,
                    "http://www.broadinstitute.org/igv/sequence/hg18", false);
        }
        return DEFAULT;
    }

    public GenomeListItem getDefaultGenomeListItem() {
        GenomeDescriptor desc = getDefaultGenomeDescriptor();
        return new GenomeListItem(desc.getName(), null, desc.getId(), 0, false);
    }

    public List<String> getChromosomeNames() {
        return genome == null ? new ArrayList() : genome.getChromosomeNames();
    }

    /**
     * A container for specific genome information which can be used to
     * manage loaded genomes.
     */
    public static class GenomeListItem {

        private int version;
        private String displayableName;
        private String location;
        private String id;
        private boolean isUserDefined = false;

        /**
         * Constructor.
         *
         * @param displayableName The name that can be shown to a user.
         * @param url             The url of the genome archive.
         * @param id              The id of the genome.
         * @param isUserDefined
         */
        public GenomeListItem(String displayableName, String url, String id, int version, boolean isUserDefined) {

            this.displayableName = displayableName;
            this.location = url;
            this.id = id;
            this.version = version;
            this.isUserDefined = isUserDefined;
        }

        public String getDisplayableName() {
            return displayableName;
        }


        public String getId() {
            return id;
        }


        public String getLocation() {
            return location;
        }


        public boolean isUserDefined() {
            return isUserDefined;
        }


        @Override
        public String toString() {
            return getDisplayableName();
        }

        /**
         * Method description
         *
         * @return
         */
        @Override
        public int hashCode() {

            int hash = 1;
            hash = hash * 31 + ((displayableName == null) ? 0 : displayableName.trim().hashCode());
            hash = hash * 13 + ((id == null) ? 0 : id.trim().hashCode());
            return hash;
        }

        /**
         * Equals method.  Two GenomeListItems are equal if their ids are equal
         *
         * @param object
         * @return
         */
        @Override
        public boolean equals(Object object) {

            if (!(object instanceof GenomeListItem)) {
                return false;
            }

            GenomeListItem item = (GenomeListItem) object;

            return getId().equals(item.getId());
        }

    }


    //////// PORTED FROM ReferenceFrame

    // /////////////////////////////////////////////////////////////////////////
    // Genome methods /////////////////////////////////////////////////////////

    /**
     * Attempt to switch genomes to newGenome.  Return the actual genome, which
     * might differ if the load of newGenome fails.
     *
     * @param newGenome
     */
    public String setGenomeId(String newGenome) {

        if (genomeId != null && !genomeId.equals(newGenome)) {
            IGVMainFrame.getInstance().getSession().getHistory().clear();
        }

        boolean loadFailed = false;

        genomeId = newGenome;
        if (!isGenomeLoaded(genomeId)) {
            try {
                if (log.isDebugEnabled()) {
                    log.debug("findGenomeAndLoad: " + genomeId);
                }
                findGenomeAndLoad(genomeId);
            } catch (IOException e) {
                log.error("Error loading genome: " + genomeId, e);
                loadFailed = true;
                MessageUtils.showMessage("Load of genome: " + genomeId + " failed.");
            }
        }
        genome = getGenome(genomeId);

        if (genome == null || loadFailed) {
            GenomeDescriptor defaultDesc = getDefaultGenomeDescriptor();
            String msg = "Could not locate genome: " + genomeId + ".  Loading " + defaultDesc.getName();
            MessageUtils.showMessage(msg);
            log.error("Could not locate genome: " + genomeId + ".  Loading " + defaultDesc.getName());

            // The previously used genome is unavailable, load the default genome
            genomeId = defaultDesc.getId();
            try {
                findGenomeAndLoad(genomeId);
            } catch (IOException e) {
                log.error("Error loading genome: " + genomeId, e);
                MessageUtils.showMessage("<html>Load of genome: " + genomeId + " failed." +
                        "<br>IGV is in an ustable state and will be closed." +
                        "<br>Please report this error to igv-help@broadinstitute.org");
                System.exit(-1);
            }

            genome = getGenome(genomeId);

            // Make this the default genome (genome loaded on startup)
            PreferenceManager.getInstance().setDefaultGenome(genomeId);
        }


        // Reset the frame manager
        IGVMainFrame.getInstance().chromosomeChangeEvent();
        return genomeId;
    }

    public String getGenomeId() {
        return genomeId;
    }

    public Genome getGenome() {
        return genome;
    }

    /**
     * Return the "home chromosome" for the current genome.
     *
     * @return
     */
    public String getHomeChr() {
        return genome == null ? Globals.CHR_ALL : genome.getHomeChromosome();
    }


}
