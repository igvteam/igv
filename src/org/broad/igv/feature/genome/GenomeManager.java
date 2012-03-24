/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.CytoBandFileParser;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.Utilities;

import java.io.*;
import java.net.MalformedURLException;
import java.net.SocketException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;

/**
 * @author jrobinso
 */
public class GenomeManager {

    private static Logger log = Logger.getLogger(GenomeManager.class);

    final public static String USER_DEFINED_GENOME_LIST_FILE = "user-defined-genomes.txt";

    private Genome currentGenome;

    private List<GenomeListItem> userDefinedGenomeArchiveList;
    private List<GenomeListItem> cachedGenomeArchiveList;
    private List<GenomeListItem> serverGenomeArchiveList;

    /**
     * The IGV instance that owns this GenomeManager.  Can be null.
     */
    IGV igv;

    public GenomeManager(IGV igv) {
        // genomeDescriptorMap = new HashMap();
        this.igv = igv;
    }

    public GenomeManager() {
        //  genomeDescriptorMap = new HashMap();
        this.igv = null;
    }

    public void setCurrentGenome(Genome currentGenome) {
        this.currentGenome = currentGenome;
    }

    /**
     * Load a genome from the given path.  Could be a .genome, or fasta file
     *
     * @param genomePath File, http, or ftp path to the .genome or indexed fasta file
     * @param monitor    ProgressMonitor  Monitor object, can be null
     * @return Genome
     * @throws FileNotFoundException
     */
    public Genome loadGenome(
            String genomePath,
            ProgressMonitor monitor)
            throws IOException {

        try {
            log.info("Loading genome: " + genomePath);
            GenomeDescriptor genomeDescriptor = null;

            if (monitor != null) {
                monitor.fireProgressChange(25);
            }

            if (genomePath.endsWith(".genome")) {

                File archiveFile = getArchiveFile(genomePath);

                genomeDescriptor = parseGenomeArchiveFile(archiveFile);
                LinkedHashMap<String, Chromosome> chromMap = loadCytobandFile(genomeDescriptor);
                Map<String, String> aliases = loadAliasFile(genomeDescriptor);

                final String id = genomeDescriptor.getId();
                final String displayName = genomeDescriptor.getName();
                boolean isFasta = genomeDescriptor.isFasta();
                currentGenome = new GenomeImpl(id, displayName, genomeDescriptor.getSequenceLocation(), isFasta);
                log.info("Genome loaded.  id= " + id);
                currentGenome.setChromosomeMap(chromMap, genomeDescriptor.isChromosomesAreOrdered());
                if (aliases != null) currentGenome.addChrAliases(aliases);
                if (!Globals.isHeadless()) {
                    updateGeneTrack(genomeDescriptor);
                }


            } else if (genomePath.endsWith(Globals.FASTA_GZIP_FILE_EXTENSION)) {
                //Assume we are dealing with a gzipped file
                boolean unzip = true;
                if (Globals.isHeadless()) {
                    log.info("Ungzipping file " + genomePath);
                } else {
                    unzip = MessageUtils.confirm("IGV cannot read gzipped fasta files directly. Would you like to un-gzip this file?");
                }

                if (!unzip) {
                    if (monitor != null) {
                        monitor.fireProgressChange(100);
                    }
                    log.info("Loading cancelled: " + genomePath);
                    return null;
                } else {
                    //This ensures we have a local copy
                    File archiveFile = getArchiveFile(genomePath);
                    String outfile = CompressionUtils.ungzipFile(archiveFile.getAbsolutePath(), null);
                    monitor.fireProgressChange(-25);
                    return loadGenome(outfile, monitor);
                }

            } else {

                // Assume its a fasta
                String fastaPath = null;
                String fastaIndexPath = null;
                if (genomePath.endsWith(".fai")) {
                    fastaPath = genomePath.substring(0, genomePath.length() - 4);
                    fastaIndexPath = genomePath;
                } else {
                    fastaPath = genomePath;
                    fastaIndexPath = genomePath + ".fai";
                }

                if (!FileUtils.resourceExists(fastaIndexPath)) {
                    //Have to make sure we have a local copy of the fasta file
                    //to index it
                    File archiveFile = getArchiveFile(fastaPath);
                    fastaPath = archiveFile.getAbsolutePath();
                    fastaIndexPath = fastaPath + ".fai";

                    log.info("Creating index file at " + fastaIndexPath);
                    FastaSequenceIndex.createIndexFile(fastaPath, fastaIndexPath);

                }

                String id = fastaPath;
                String name = (new File(fastaPath)).getName();
                if (HttpUtils.isURL(fastaPath)) {
                    name = Utilities.getFileNameFromURL(fastaPath);
                }

                currentGenome = new GenomeImpl(id, name, fastaPath, true);

                log.info("Genome loaded.  id= " + id);
                if (!Globals.isHeadless()) {
                    IGV.getInstance().createGeneTrack(currentGenome, null, null, null, null);
                }
            }

            if (monitor != null) {
                monitor.fireProgressChange(25);
            }

            // Do this last so that user defined aliases have preference.
            currentGenome.loadUserDefinedAliases();

            return currentGenome;

        } catch (SocketException e) {
            throw new GenomeServerException("Server connection error", e);
        }

    }

    /**
     * Returns a File of the provided genomePath.
     * If the genomePath is a URL, it will be downloaded
     * and saved in the genome cache directory.
     *
     * @param genomePath
     * @return
     * @throws MalformedURLException
     * @throws UnsupportedEncodingException
     */
    private File getArchiveFile(String genomePath) throws MalformedURLException, UnsupportedEncodingException {
        File archiveFile;
        if (HttpUtils.getInstance().isURL(genomePath.toLowerCase())) {
            // We need a local copy, as there is no http zip file reader
            URL genomeArchiveURL = new URL(genomePath);
            String cachedFilename = Utilities.getFileNameFromURL(
                    URLDecoder.decode(new URL(genomePath).getFile(), "UTF-8"));
            if (!DirectoryManager.getGenomeCacheDirectory().exists()) {
                DirectoryManager.getGenomeCacheDirectory().mkdir();
            }
            archiveFile = new File(DirectoryManager.getGenomeCacheDirectory(), cachedFilename);
            refreshCache(archiveFile, genomeArchiveURL);
        } else {
            archiveFile = new File(genomePath);
        }
        return archiveFile;
    }

    private void updateGeneTrack(GenomeDescriptor genomeDescriptor) throws IOException {

        InputStream geneStream = null;
        try {
            geneStream = genomeDescriptor.getGeneStream();
            BufferedReader reader = geneStream == null ? null : new BufferedReader(new InputStreamReader(geneStream));
            IGV.getInstance().createGeneTrack(currentGenome, reader,
                    genomeDescriptor.getGeneFileName(), genomeDescriptor.getGeneTrackName(),
                    genomeDescriptor.getUrl());
        } finally {
            if (geneStream != null) geneStream.close();
        }
    }

    /**
     * Load the cytoband file specified in the genome descriptor and return an ordered hash map of
     * chromsome name -> chromosome.   This is a legacy method, kept for backward compatibiltiy of
     * .genome files in which the chromosome lengths are specified in a cytoband file.
     *
     * @param genomeDescriptor
     * @return
     */
    private LinkedHashMap<String, Chromosome> loadCytobandFile(GenomeDescriptor genomeDescriptor) {
        InputStream is = null;
        try {

            is = genomeDescriptor.getCytoBandStream();
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));
            return CytoBandFileParser.loadData(reader);

        } catch (IOException ex) {
            log.warn("Error loading cytoband file", ex);
            throw new RuntimeException("Error loading cytoband file" + genomeDescriptor.cytoBandFileName);
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

    /**
     * Load the chromosome alias file, if any, specified in the genome descriptor.
     *
     * @param genomeDescriptor
     * @return The chromosome alias map, or null if none is defined.
     */
    private Map<String, String> loadAliasFile(GenomeDescriptor genomeDescriptor) {
        InputStream aliasStream = null;
        try {
            aliasStream = genomeDescriptor.getChrAliasStream();
            if (aliasStream != null) {
                Map<String, String> chrAliasTable = new HashMap();
                BufferedReader reader = new BufferedReader(new InputStreamReader(aliasStream));
                String nextLine = "";
                while ((nextLine = reader.readLine()) != null) {
                    String[] kv = nextLine.split("\t");
                    if (kv.length > 1) {
                        chrAliasTable.put(kv[0], kv[1]);
                    }
                }

                return chrAliasTable;
            } else {
                return null;
            }
        } catch (Exception e) {
            // We don't want to bomb if the alias load fails.  Just log it and proceed.
            log.error("Error loading chromosome alias table");
            return null;
        } finally {
            try {
                if (aliasStream != null) {
                    aliasStream.close();
                }
            } catch (IOException ex) {
                log.warn("Error closing zip stream!", ex);
            }
        }
    }

    /**
     * Refresh a locally cached genome
     *
     * @param cachedFile
     * @param genomeArchiveURL
     * @throws IOException
     */
    private void refreshCache(File cachedFile, URL genomeArchiveURL) {
        // Look in cache first


        try {
            if (cachedFile.exists()) {


                boolean remoteModfied = !HttpUtils.getInstance().compareResources(cachedFile, genomeArchiveURL);

                // Force an update of cached genome if file length does not equal remote content length
                boolean forceUpdate = remoteModfied &&
                        PreferenceManager.getInstance().getAsBoolean(PreferenceManager.AUTO_UPDATE_GENOMES);
                if (forceUpdate) {
                    log.info("Refreshing genome: " + genomeArchiveURL.toString());
                    File tmpFile = new File(cachedFile.getAbsolutePath() + ".tmp");
                    if (HttpUtils.getInstance().downloadFile(genomeArchiveURL.toExternalForm(), tmpFile)) {
                        FileUtils.copyFile(tmpFile, cachedFile);
                        tmpFile.deleteOnExit();
                    }
                }

            } else {
                // Copy file directly from the server to local cache.
                HttpUtils.getInstance().downloadFile(genomeArchiveURL.toExternalForm(), cachedFile);
            }
        } catch (Exception e) {
            log.error("Error refreshing genome cache. ", e);
            MessageUtils.showMessage(("An error was encountered refreshing the genome cache: " + e.getMessage() +
                    "<br> If this problem persists please contact igv-help@broadinstitute.org"));
        }

    }


    /**
     * Creates a genome descriptor.
     */
    public GenomeDescriptor parseGenomeArchiveFile(File f)
            throws IOException {


        String zipFilePath = f.getAbsolutePath();

        if (!f.exists()) {
            log.error("Genome file: " + f.getAbsolutePath() + " does not exist.");
            return null;
        }


        GenomeDescriptor genomeDescriptor = null;
        Map<String, ZipEntry> zipEntries = new HashMap();
        ZipFile zipFile = new ZipFile(zipFilePath);

        FileInputStream fileInputStream = null;
        try {
            fileInputStream = new FileInputStream(f);
            ZipInputStream zipInputStream = new ZipInputStream(fileInputStream);
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

                    if ((sequenceLocation != null) && !HttpUtils.getInstance().isURL(sequenceLocation)) {
                        File sequenceFolder = null;
                        // Relative or absolute location?
                        if (sequenceLocation.startsWith("/") || sequenceLocation.startsWith("\\")) {
                            sequenceFolder = new File(sequenceLocation);
                        } else {
                            File tempZipFile = new File(zipFilePath);
                            sequenceFolder = new File(tempZipFile.getParent(), sequenceLocation);

                        }
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

                    boolean fasta = false;
                    String fastaString = properties.getProperty("fasta");
                    if (fastaString != null) {
                        try {
                            fasta = Boolean.parseBoolean(fastaString);
                        } catch (Exception e) {
                            log.error("Error parsing version string: " + versionString);
                        }
                    }

                    boolean chromosomesAreOrdered = false;
                    String tmp = properties.getProperty(Globals.GENOME_ORDERED_KEY);
                    if (tmp != null) {
                        try {
                            chromosomesAreOrdered = Boolean.parseBoolean(tmp);
                        } catch (Exception e) {
                            log.error("Error parsing ordered string: " + tmp);
                        }
                    }

                    String url = properties.getProperty(Globals.GENOME_URL_KEY);


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
                            chromosomesAreOrdered,
                            fasta);

                    if (url != null) {
                        genomeDescriptor.setUrl(url);
                    }

                }
                zipEntry = zipInputStream.getNextEntry();
            }
        } finally {
            try {
                if (fileInputStream != null) {
                    fileInputStream.close();
                }
            } catch (IOException ex) {
                log.warn("Error closing imported genome zip stream!", ex);
            }
        }
        return genomeDescriptor;
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
    public List<GenomeListItem> getServerGenomeArchiveList(Set excludedArchivesUrls)
            throws IOException {

        if (serverGenomeListUnreachable) {
            return null;
        }

        if (serverGenomeArchiveList == null) {
            serverGenomeArchiveList = new LinkedList();
            BufferedReader dataReader = null;
            InputStream inputStream = null;
            String genomeListURLString = "";
            try {
                genomeListURLString = PreferenceManager.getInstance().getGenomeListURL();
                URL serverGenomeURL = new URL(genomeListURLString);

                if (HttpUtils.getInstance().isURL(genomeListURLString)) {
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

                            //String displayableName, String url, String id, int version, boolean isUserDefined
                            String name = fields[0];
                            String url = fields[1];
                            String id = fields[2];

                            boolean valid = true;
                            if (url.length() == 0) {
                                log.error("Genome entry : " + name + " has an empty URL string.  Check for extra tabs in the definition file: " +
                                        PreferenceManager.getInstance().getGenomeListURL());
                                valid = false;
                            }
                            // TODO -- more validation


                            if (valid) {


                                try {
                                    GenomeListItem item = new GenomeListItem(fields[0], fields[1], fields[2], false);
                                    serverGenomeArchiveList.add(item);
                                } catch (Exception e) {
                                    log.error(
                                            "Error reading a line from server genome list" + " line was: [" + genomeRecord + "]",
                                            e);
                                }
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
                        genomeListURLString + ").    Only locally defined genomes will be available.",
                        PreferenceManager.SHOW_GENOME_SERVER_WARNING);
            } finally {
                if (dataReader != null) {
                    dataReader.close();
                }
                if (inputStream != null) {
                    inputStream.close();
                }
            }
        }
        return serverGenomeArchiveList;
    }

    /**
     * Gets a list of all the user-defined genome archive files that
     * IGV knows about.
     *
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public List<GenomeListItem> getUserDefinedGenomeArchiveList()
            throws IOException {


        if (userDefinedGenomeArchiveList == null) {

            boolean updateClientGenomeListFile = false;

            userDefinedGenomeArchiveList = new LinkedList();

            File listFile = new File(DirectoryManager.getGenomeCacheDirectory(), USER_DEFINED_GENOME_LIST_FILE);

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

                    GenomeListItem item = new GenomeListItem(fields[0], file, fields[2], true);
                    userDefinedGenomeArchiveList.add(item);
                }
            } finally {
                if (reader != null) reader.close();
            }
            if (updateClientGenomeListFile) {
                updateImportedGenomePropertyFile();
            }
        }
        return userDefinedGenomeArchiveList;
    }

    /**
     * Method description
     */
    public void clearGenomeCache() {

        File[] files = DirectoryManager.getGenomeCacheDirectory().listFiles();
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
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public List<GenomeListItem> getCachedGenomeArchiveList()
            throws IOException {

        if (cachedGenomeArchiveList == null) {
            cachedGenomeArchiveList = new LinkedList();

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
                                    false);
                    cachedGenomeArchiveList.add(item);
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
        }

        return cachedGenomeArchiveList;
    }


    /**
     * Reconstructs the user-define genome property file.
     *
     * @throws IOException
     */
    public void updateImportedGenomePropertyFile() {

        if (userDefinedGenomeArchiveList == null) {
            return;
        }

        File listFile = new File(DirectoryManager.getGenomeCacheDirectory(), USER_DEFINED_GENOME_LIST_FILE);
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
            for (GenomeListItem genomeListItem : userDefinedGenomeArchiveList) {
                writer.print(genomeListItem.getDisplayableName());
                writer.print("\t");
                writer.print(genomeListItem.getLocation());
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
            log.error("Error updating genome property file", e);
            MessageUtils.showMessage("Error updating user-defined genome list " + e.getMessage());

        } finally {
            if (writer != null) writer.close();
            if (backup != null) backup.delete();
        }
    }

    /**
     * Create a genome archive (.genome) file.
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

        if (monitor != null) monitor.fireProgressChange(25);

        File archiveFile = (new GenomeImporter()).createGenomeArchive(zipFileLocation,
                genomeFileName, genomeId, genomeDisplayName, relativeSequenceLocation,
                fastaInputFile, refFlatFile, cytobandFile, chrAliasFile,
                sequenceOutputLocationOverride, monitor);

        if (monitor != null) monitor.fireProgressChange(75);

        if (archiveFile == null) {
            return null;
        } else {
            GenomeListItem newItem = new GenomeListItem(genomeDisplayName, archiveFile.getAbsolutePath(), genomeId, true);
            addUserDefineGenomeItem(newItem);
            return newItem;
        }
    }

    public String getGenomeId() {
        return currentGenome == null ? null : currentGenome.getId();
    }

    public Genome getCurrentGenome() {
        return currentGenome;
    }

    public void addUserDefineGenomeItem(GenomeListItem genomeListItem) {
        userDefinedGenomeArchiveList.add(0, genomeListItem);
        updateImportedGenomePropertyFile();
    }
}
