/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
import org.broad.igv.dev.affective.AffectiveUtils;
import org.broad.igv.feature.*;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.broad.igv.util.collections.CI;

import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.SocketException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.*;
import java.util.List;
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

    private static GenomeManager theInstance;

    private Genome currentGenome;

    private LinkedHashSet<GenomeListItem> userDefinedGenomeArchiveList;
    private List<GenomeListItem> serverGenomeArchiveList;
    private List<GenomeListItem> cachedGenomeArchiveList;
    private Set<String> excludedArchivesUrls = new HashSet();

    //ID comparison will be case insensitive
    private Map<String, GenomeListItem> genomeItemMap = new CI.CILinkedHashMap<GenomeListItem>();

    public static void main(String[] args) {
        if (args.length >= 1 && args[0].equals("genList")) {
            if (args.length != 4)
                throw new IllegalArgumentException("Incorrect number of inputs, expected genList [dir] [rootPath] [outFile]");
        }
    }


    public synchronized static GenomeManager getInstance() {
        if (theInstance == null) {
            theInstance = new GenomeManager();
        }
        return theInstance;
    }

    private GenomeManager() {

    }


    public void setCurrentGenome(Genome currentGenome) {
        if (currentGenome != null) {
            PreferenceManager.getInstance().setDefaultGenome(currentGenome.getId());
        }
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

            GenomeImpl newGenome = null;

            if (monitor != null) {
                monitor.fireProgressChange(25);
            }

            if (genomePath.endsWith(Globals.GZIP_FILE_EXTENSION)) {
                throw new GenomeException("IGV cannot readed gzipped genome files.  Please un-gzip the file and try again.");
            } else if (genomePath.endsWith(".genome")) {
                newGenome = loadDotGenomeFile(genomePath);
            } else if (genomePath.endsWith(".gbk")) {
                newGenome = loadGenbankFile(genomePath);
            } else {
                // Assume a fasta file
                newGenome = loadFastaFile(genomePath);
            }

            if (monitor != null) {
                monitor.fireProgressChange(25);
            }

            setCurrentGenome(newGenome);

            if (IGV.hasInstance() && !Globals.isHeadless()) {
                FeatureTrack geneFeatureTrack = newGenome.getGeneTrack();
                IGV.getInstance().setGenomeTracks(geneFeatureTrack);
            }

            log.info("Genome loaded.  id= " + newGenome.getId());

            return currentGenome;

        } catch (SocketException e) {
            throw new GenomeServerException("Server connection error", e);
        }

    }

    private GenomeImpl loadGenbankFile(String genomePath) throws IOException {
        GenomeImpl newGenome;
        GenbankParser genbankParser = new GenbankParser(genomePath);

        String chr = genbankParser.getAccession();
        String name = genbankParser.getLocusName();
        if (!name.equals(chr)) {
            name = name + " (" + chr + ")";
        }

        byte[] seq = genbankParser.getSequence();
        Sequence sequence = new InMemorySequence(chr, seq);
        newGenome = new GenomeImpl(chr, name, sequence);
        newGenome.loadUserDefinedAliases();
        setCurrentGenome(newGenome);

        if (IGV.hasInstance() && !Globals.isHeadless()) {
            FeatureTrack geneFeatureTrack = createGeneTrack(newGenome, genbankParser.getFeatures());
            newGenome.setGeneTrack(geneFeatureTrack);
        }

        return newGenome;
    }

    /**
     * Create a Genome from a single fasta file.
     *
     * @param genomePath
     * @return
     * @throws IOException
     */
    private GenomeImpl loadFastaFile(String genomePath) throws IOException {
        GenomeImpl newGenome;// Assume its a fasta
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
            if (!FileUtils.isRemote(fastaPath)) {
                File archiveFile = getArchiveFile(fastaPath);
                fastaPath = archiveFile.getAbsolutePath();
                fastaIndexPath = fastaPath + ".fai";

                FastaUtils.createIndexFile(fastaPath, fastaIndexPath);
            }

        }

        GenomeListItem item = buildFromPath(fastaPath);
        if (item == null) {
            throw new IOException(fastaPath + " does not exist, could not load genome");
        }

        FastaIndexedSequence fastaSequence = new FastaIndexedSequence(fastaPath);
        Sequence sequence = new SequenceWrapper(fastaSequence);
        newGenome = new GenomeImpl(item.getId(), item.getDisplayableName(), sequence);
        newGenome.loadUserDefinedAliases();
        setCurrentGenome(newGenome);
        return newGenome;
    }

    /**
     * @param path
     * @return GenomeListItem representing this path, or null
     *         if it's a local path which doesn't exist (we don't check for existence of remote file)
     */
    public static GenomeListItem buildFromPath(String path) {

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
     * Create a genome from a ".genome" file.  In addition to the reference sequence .genome files can optionally
     * specify cytobands and annotations.
     *
     * @param genomePath
     * @return
     * @throws IOException
     */
    private GenomeImpl loadDotGenomeFile(String genomePath) throws IOException {
        GenomeImpl newGenome;
        File archiveFile = getArchiveFile(genomePath);

        GenomeDescriptor genomeDescriptor = parseGenomeArchiveFile(archiveFile);

        final String id = genomeDescriptor.getId();
        final String displayName = genomeDescriptor.getName();

        boolean isFasta = genomeDescriptor.isFasta();
        String[] fastaFiles = genomeDescriptor.getFastaFileNames();

        LinkedHashMap<String, List<Cytoband>> cytobandMap = null;
        if (genomeDescriptor.hasCytobands()) {
            cytobandMap = loadCytobandFile(genomeDescriptor);
        }


        String sequencePath = genomeDescriptor.getSequenceLocation();
        Sequence sequence = null;
        //We preserve ordering only for legacy genomes
        boolean chromosOrdered = false;
        if (sequencePath == null) {
            sequence = null;
        } else if (!isFasta) {
            // Legacy genomes
            sequencePath = SequenceWrapper.checkSequenceURL(sequencePath);
            IGVSequence igvSequence = new IGVSequence(sequencePath);
            if (cytobandMap != null) {
                chromosOrdered = genomeDescriptor.isChromosomesAreOrdered();
                igvSequence.generateChromosomes(cytobandMap, chromosOrdered);
            }
            sequence = new SequenceWrapper(igvSequence);
        } else if (fastaFiles != null) {
            FastaDirectorySequence fastaDirectorySequence = new FastaDirectorySequence(sequencePath, fastaFiles);
            sequence = new SequenceWrapper(fastaDirectorySequence);
        } else {
            FastaIndexedSequence fastaSequence = new FastaIndexedSequence(sequencePath);
            sequence = new SequenceWrapper(fastaSequence);
        }

        newGenome = new GenomeImpl(id, displayName, sequence, chromosOrdered);
        if (cytobandMap != null) {
            newGenome.setCytobands(cytobandMap);
        }

        Collection<Collection<String>> aliases = loadChrAliases(genomeDescriptor);
        if (aliases != null) {
            newGenome.addChrAliases(aliases);
        }
        // Do this last so that user defined aliases have preference.
        newGenome.loadUserDefinedAliases();
        setCurrentGenome(newGenome);


        InputStream geneStream = null;
        if (genomeDescriptor.getGeneFileName() != null) {
            try {
                geneStream = genomeDescriptor.getGeneStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(geneStream));
                FeatureTrack geneFeatureTrack = createGeneTrack(newGenome, reader,
                        genomeDescriptor.getGeneFileName(), genomeDescriptor.getGeneTrackName(),
                        genomeDescriptor.getUrl());

                newGenome.setGeneTrack(geneFeatureTrack);
            } finally {
                if (geneStream != null) geneStream.close();
            }
        }

        genomeDescriptor.close();
        return newGenome;
    }

    /**
     * Returns a File of the provided genomePath. If the genomePath is a URL, it will be downloaded
     * and saved in the genome cache directory.
     *
     * @param genomePath
     * @return
     * @throws MalformedURLException
     * @throws UnsupportedEncodingException
     */
    private File getArchiveFile(String genomePath) throws MalformedURLException, UnsupportedEncodingException {
        File archiveFile;
        if (HttpUtils.isRemoteURL(genomePath.toLowerCase())) {
            // We need a local copy, as there is no http zip file reader
            URL genomeArchiveURL = new URL(genomePath);
            final String tmp = URLDecoder.decode(new URL(genomePath).getFile(), "UTF-8");
            String cachedFilename = Utilities.getFileNameFromURL(tmp);
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


    /**
     * Load the cytoband file specified in the genome descriptor.
     *
     * @param genomeDescriptor
     * @return Cytobands as a map keyed by chromosome
     */
    private LinkedHashMap<String, List<Cytoband>> loadCytobandFile(GenomeDescriptor genomeDescriptor) {
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

    static Collection<Collection<String>> loadChrAliases(BufferedReader br) throws IOException {
        String nextLine = "";
        Collection<Collection<String>> synonymList = new ArrayList<Collection<String>>();
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            if (tokens.length > 1) {
                Collection<String> synonyms = new ArrayList<String>();
                for (String t : tokens) {
                    String syn = t.trim();
                    if (t.length() > 0) synonyms.add(syn.trim());
                }
                synonymList.add(synonyms);
            }
        }
        return synonymList;
    }

    /**
     * Load the chromosome alias file, if any, specified in the genome descriptor.
     *
     * @param genomeDescriptor
     * @return The chromosome alias map, or null if none is defined.
     */
    private Collection<Collection<String>> loadChrAliases(GenomeDescriptor genomeDescriptor) {
        InputStream aliasStream = null;
        try {
            aliasStream = genomeDescriptor.getChrAliasStream();
            if (aliasStream != null) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(aliasStream));
                return loadChrAliases(reader);
            } else {
                return null;
            }
        } catch (IOException e) {
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
            MessageUtils.showErrorMessage("An error was encountered refreshing the genome cache: " + e.getMessage(), e);
        }

    }


    /**
     * Creates a genome descriptor.
     */
    public GenomeDescriptor parseGenomeArchiveFile(File f)
            throws IOException {


        if (!f.exists()) {
            throw new FileNotFoundException("Genome file: " + f.getAbsolutePath() + " does not exist.");
        }

        GenomeDescriptor genomeDescriptor = null;
        Map<String, ZipEntry> zipEntries = new HashMap();
        ZipFile zipFile = new ZipFile(f);

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

                    String cytobandZipEntryName = properties.getProperty(Globals.GENOME_ARCHIVE_CYTOBAND_FILE_KEY);
                    String geneFileName = properties.getProperty(Globals.GENOME_ARCHIVE_GENE_FILE_KEY);
                    String chrAliasFileName = properties.getProperty(Globals.GENOME_CHR_ALIAS_FILE_KEY);
                    String sequenceLocation = properties.getProperty(Globals.GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY);

                    if ((sequenceLocation != null) && !HttpUtils.isRemoteURL(sequenceLocation)) {
                        File sequenceFolder = null;
                        // Relative or absolute location? We use a few redundant methods to check,
                        //since we don't know what platform the file was created on or is running on
                        sequenceFolder = new File(sequenceLocation);
                        boolean isAbsolutePath = sequenceFolder.isAbsolute() ||
                                sequenceLocation.startsWith("/") || sequenceLocation.startsWith("\\");
                        if (!isAbsolutePath) {
                            sequenceFolder = new File(f.getParent(), sequenceLocation);
                        }
                        sequenceLocation = sequenceFolder.getCanonicalPath();
                        sequenceLocation.replace('\\', '/');
                    }

                    boolean chrNamesAltered = parseBooleanPropertySafe(properties, "filenamesAltered");
                    boolean fasta = parseBooleanPropertySafe(properties, "fasta");
                    boolean fastaDirectory = parseBooleanPropertySafe(properties, "fastaDirectory");
                    boolean chromosomesAreOrdered = parseBooleanPropertySafe(properties, Globals.GENOME_ORDERED_KEY);


                    String fastaFileNameString = properties.getProperty("fastaFiles");
                    String url = properties.getProperty(Globals.GENOME_URL_KEY);


                    // The new descriptor
                    genomeDescriptor = new GenomeZipDescriptor(
                            properties.getProperty(Globals.GENOME_ARCHIVE_NAME_KEY),
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
                            fasta,
                            fastaDirectory,
                            fastaFileNameString);

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

    private boolean parseBooleanPropertySafe(Properties properties, String key) {
        String propertyString = properties.getProperty(key);
        return Boolean.parseBoolean(propertyString);
    }

    boolean serverGenomeListUnreachable = false;

    /**
     * Get the server genome list, if accessible. If not,
     * get the list of cached genomes
     *
     * @return
     * @throws IOException
     */
    public List<GenomeListItem> getGenomeArchiveList() {
        List<GenomeListItem> genomeArchiveList = getServerGenomeArchiveList(this.excludedArchivesUrls);

        if (genomeArchiveList == null) {
            try {
                genomeArchiveList = getCachedGenomeArchiveList();
            } catch (IOException e) {
                MessageUtils.showErrorMessage("Cannot access cached genome list", e);
            }
        }
        return genomeArchiveList;
    }

    /**
     * Check the server or cache for the given {@code genomeID}, load it into the current set.
     *
     * @param genomeId
     * @return True if found, false if not
     * @throws IOException
     */
    public boolean loadFromArchive(String genomeId) throws IOException {
        GenomeListItem matchingItem = findGenomeListItemById(genomeId);
        if (matchingItem != null) {
            GenomeManager.getInstance().addGenomeItems(Arrays.asList(matchingItem));
        }
        return matchingItem != null;
    }

    /**
     * Gets a list of all the server genome archive files that
     * IGV knows about.
     *
     * @param excludedArchivesUrls The set of file location to exclude in the return list.
     * @return List<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public List<GenomeListItem> getServerGenomeArchiveList(Set excludedArchivesUrls) {

        if (serverGenomeListUnreachable) {
            return null;
        }

        if (serverGenomeArchiveList == null) {
            serverGenomeArchiveList = new LinkedList<GenomeListItem>();
            BufferedReader dataReader = null;
            InputStream inputStream = null;
            String genomeListURLString = "";
            try {
                genomeListURLString = PreferenceManager.getInstance().getGenomeListURL();
                URL serverGenomeURL = new URL(genomeListURLString);

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

                            // Throw away records we don't want to see
                            if (excludedArchivesUrls != null) {
                                if (excludedArchivesUrls.contains(fields[1])) {
                                    continue;
                                }
                            }

                            String name = fields[0];
                            String url = fields[1];
                            String id = fields[2];
                            GenomeListItem item = new GenomeListItem(name, url, id);
                            serverGenomeArchiveList.add(item);

                        } else {
                            log.error("Found invalid server genome list record: " + genomeRecord);
                        }
                    }
                }
            } catch (Exception e) {
                serverGenomeListUnreachable = true;
                serverGenomeArchiveList = null;
                log.error("Error fetching genome list: ", e);
                ConfirmDialog.optionallyShowInfoDialog("Warning: could not connect to the genome server (" +
                        genomeListURLString + ").    Only locally defined genomes will be available.",
                        PreferenceManager.SHOW_GENOME_SERVER_WARNING);
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

        return serverGenomeArchiveList;
    }

    /**
     * Searches through currently loaded GenomeListItems and returns
     * that with a matching ID. null if not found. To search through
     * all server and user defined genomes, use #findGenomeListItemByItem
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
    public GenomeListItem findGenomeListItemById(String genomeId) throws IOException {
        GenomeListItem matchingItem = genomeItemMap.get(genomeId);
        if (matchingItem == null) {
            // If genome archive was not found, check things not
            //currently loaded
            matchingItem = searchGenomeList(genomeId, getGenomeArchiveList());
            if (matchingItem != null) return matchingItem;

            matchingItem = searchGenomeList(genomeId, getUserDefinedGenomeArchiveList());
            if (matchingItem != null) return matchingItem;

        }
        return matchingItem;
    }

    private GenomeListItem searchGenomeList(String genomeId, Iterable<GenomeListItem> genomeList) {
        for (GenomeListItem item : genomeList) {
            if (item.getId().equals(genomeId)) {
                return item;
            }
        }
        return null;
    }

    public Collection<GenomeListItem> getGenomes() {
        return genomeItemMap.values();
    }

    /**
     * Completely rebuild the genome drop down info.
     * This will load genomes from server/cached/user defined (as appropriate),
     * and only keep the ones that the user has chosen to be displayed
     */
    public void buildGenomeItemList() {

        // Build a single available genome list from both client, server
        // and cached information. This allows us to process
        // everything the same way.
        Collection<GenomeListItem> tmpuserDefinedGenomeList = null;
        Collection<GenomeListItem> tmpArchiveGenomeItemList = null;

        boolean affectiveMode = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.AFFECTIVE_ENABLE);
        if (affectiveMode) {
            tmpArchiveGenomeItemList = Arrays.asList(AffectiveUtils.GENOME_DESCRIPTOR);
        } else {
            tmpArchiveGenomeItemList = getGenomeArchiveList();
            try {
                tmpuserDefinedGenomeList = getUserDefinedGenomeArchiveList();
            } catch (IOException e) {
                MessageUtils.showErrorMessage("Cannot access user defined genome archive list", e);
            }
        }

        combineGenomeLists(tmpuserDefinedGenomeList, tmpArchiveGenomeItemList);

    }


    /**
     * Combine our different lists of genomeListItems
     *
     * @param userDefinedGenomeItemList
     * @param archiveGenomeItemList     Either server or cached genomes
     */
    private void combineGenomeLists(Collection<GenomeListItem> userDefinedGenomeItemList,
                                    Collection<GenomeListItem> archiveGenomeItemList) {

        //We use a LinkedHashMap to prevent loading duplicates
        genomeItemMap = new LinkedHashMap<String, GenomeListItem>();

        String[] genomeIdArray = PreferenceManager.getInstance().getGenomeIdDisplayList();

        if (genomeIdArray.length == 0) {
            genomeIdArray = new String[]{PreferenceManager.getInstance().getDefaultGenome(), "hg18"};
        }

        if ((userDefinedGenomeItemList != null) && !userDefinedGenomeItemList.isEmpty()) {
            addAllInSet(userDefinedGenomeItemList, genomeIdArray, genomeItemMap);
        }

        addAllInSet(archiveGenomeItemList, genomeIdArray, genomeItemMap);
    }

    /**
     * Adds the {@code genomeListItems} to the map by their ids,
     * in order of {@code keepGenomeIds},
     * iff the id is in {@code keepGenomeIds} and NOT in genomeMap
     *
     * @param genomeListItems
     * @param keepGenomeIds
     * @param genomeMap
     */
    private void addAllInSet(Iterable<GenomeListItem> genomeListItems, String[] keepGenomeIds, Map<String, GenomeListItem> genomeMap) {
        for (String id : keepGenomeIds) {
            GenomeListItem genomeListItem = searchGenomeList(id, genomeListItems);

            //if we didn't find the id, it may be a path
            if (genomeListItem == null) {
                genomeListItem = buildFromPath(id);
            }

            if (!genomeMap.containsKey(id) && genomeListItem != null) {
                genomeMap.put(id, genomeListItem);
            }

        }
    }


    /**
     * Gets a list of all the user-defined genome archive files that
     * IGV knows about.
     *
     * @return LinkedHashSet<GenomeListItem>
     * @throws IOException
     * @see GenomeListItem
     */
    public Collection<GenomeListItem> getUserDefinedGenomeArchiveList() throws IOException {


        if (userDefinedGenomeArchiveList == null) {

            boolean updateClientGenomeListFile = false;

            userDefinedGenomeArchiveList = new LinkedHashSet<GenomeListItem>();

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

                    GenomeListItem item = new GenomeListItem(fields[0], file, fields[2]);
                    userDefinedGenomeArchiveList.add(item);
                }
            } catch (FileNotFoundException e) {
                //We swallow this because the user may not have the file,
                //which doesn't really matter
                log.info(e);
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
    private List<GenomeListItem> getCachedGenomeArchiveList()
            throws IOException {

        if (cachedGenomeArchiveList == null) {
            cachedGenomeArchiveList = new LinkedList<GenomeListItem>();

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
                                    properties.getProperty(Globals.GENOME_ARCHIVE_ID_KEY));
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
            MessageUtils.showErrorMessage("Error updating user-defined genome list " + e.getMessage(), e);

        } finally {
            if (writer != null) writer.close();
            if (backup != null) backup.delete();
        }
    }

    /**
     * Create a genome archive (.genome) file.
     *
     * @param genomeFile
     * @param cytobandFileName  A File path to a file that contains cytoband data.
     * @param refFlatFileName   A File path to a gene file.
     * @param fastaFileName     A File path to a FASTA file, a .gz file containing a
     *                          single FASTA file, or a directory containing ONLY FASTA files.
     *                          (relative to the .genome file to be created) where the sequence data for
     *                          the new genome will be written.
     * @param genomeDisplayName The unique user-readable name of the new genome.
     * @param genomeId          The id to be assigned to the genome.
     * @param monitor           A ProgressMonitor used to track progress - null,
     *                          if no progress updating is required.
     * @return GenomeListItem
     * @throws FileNotFoundException
     */
    public GenomeListItem defineGenome(File genomeFile,
                                       String cytobandFileName,
                                       String refFlatFileName,
                                       String fastaFileName,
                                       String chrAliasFileName,
                                       String genomeDisplayName,
                                       String genomeId,
                                       ProgressMonitor monitor)
            throws IOException {

        File refFlatFile = null;
        File cytobandFile = null;
        File chrAliasFile = null;

        if (genomeFile != null) {
            PreferenceManager.getInstance().setLastGenomeImportDirectory(genomeFile.getParentFile());
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

        if (monitor != null) monitor.fireProgressChange(25);

        (new GenomeImporter()).createGenomeArchive(genomeFile, genomeId,
                genomeDisplayName, fastaFileName, refFlatFile, cytobandFile, chrAliasFile);

        if (monitor != null) monitor.fireProgressChange(75);

        GenomeListItem newItem = new GenomeListItem(genomeDisplayName, genomeFile.getAbsolutePath(), genomeId);
        addGenomeItem(newItem);
        return newItem;

    }

    public String getGenomeId() {
        return currentGenome == null ? null : currentGenome.getId();
    }

    public Genome getCurrentGenome() {
        return currentGenome;
    }

    public void addGenomeItems(Collection<GenomeListItem> genomeListItems) {
        for (GenomeListItem genomeListItem : genomeListItems) {
            genomeItemMap.put(genomeListItem.getId(), genomeListItem);
            userDefinedGenomeArchiveList.add(genomeListItem);
        }
        PreferenceManager.getInstance().saveGenomeIdDisplayList(genomeItemMap.values());
        updateImportedGenomePropertyFile();
    }

    public void addGenomeItem(GenomeListItem genomeListItem) {
        genomeItemMap.put(genomeListItem.getId(), genomeListItem);
        PreferenceManager.getInstance().saveGenomeIdDisplayList(genomeItemMap.values());
        userDefinedGenomeArchiveList.add(genomeListItem);
        updateImportedGenomePropertyFile();
    }

    /**
     * Given a directory, looks for all .genome files,
     * and outputs a list of these genomes suitable for parsing by IGV.
     * Intended to be run on server periodically.
     *
     * @param inDir    Directory in which all genome files live
     * @param rootPath The path to be prepended to file names (e.g. http://igvdata.broadinstitute.org)
     * @param outPath  Path to output file, where we will write the results
     */
    public void generateGenomeList(File inDir, String rootPath, String outPath) {
        File[] genomeFiles = inDir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name == null) return false;
                return name.toLowerCase().endsWith(".genome");
            }
        });

        PrintWriter writer;
        try {
            writer = new PrintWriter(outPath);
        } catch (FileNotFoundException e) {
            log.error("Error opening " + outPath);
            e.printStackTrace();
            return;
        }


        GenomeDescriptor descriptor;
        for (File f : genomeFiles) {
            String curLine = "";
            try {
                descriptor = parseGenomeArchiveFile(f);
                curLine += descriptor.getName();
                curLine += "\t" + rootPath + "/" + f.getName();
                curLine += "\t" + descriptor.getId();
            } catch (IOException e) {
                log.error("Error parsing genome file. Skipping " + f.getAbsolutePath());
                log.error(e);
                continue;
            }
            writer.println(curLine);
        }

        writer.close();

    }

    /**
     * @param reader        a reader for the gene (annotation) file.
     * @param genome
     * @param geneFileName
     * @param geneTrackName
     */
    public FeatureTrack createGeneTrack(Genome genome, BufferedReader reader, String geneFileName, String geneTrackName,
                                        String annotationURL) {

        FeatureDB.clearFeatures();
        FeatureTrack geneFeatureTrack = null;

        if (reader != null) {
            FeatureParser parser;
            if (geneFileName.endsWith(".embl")) {
                parser = new EmblFeatureTableParser();
            } else if (GFFFeatureSource.isGFF(geneFileName)) {
                parser = new GFFParser(geneFileName);
            } else {
                parser = AbstractFeatureParser.getInstanceFor(new ResourceLocator(geneFileName), genome);
            }
            if (parser == null) {
                MessageUtils.showMessage("ERROR: Unrecognized annotation file format: " + geneFileName +
                        "<br>Annotations for genome: " + genome.getId() + " will not be loaded.");
            } else {
                List<org.broad.tribble.Feature> genes = parser.loadFeatures(reader, genome);
                String name = geneTrackName;
                if (name == null) name = "Genes";

                String id = genome.getId() + "_genes";
                geneFeatureTrack = new FeatureTrack(id, name, new FeatureCollectionSource(genes, genome));
                geneFeatureTrack.setMinimumHeight(5);
                geneFeatureTrack.setHeight(35);
                geneFeatureTrack.setTrackType(TrackType.GENE);
                geneFeatureTrack.setColor(Color.BLUE.darker());
                TrackProperties props = parser.getTrackProperties();
                if (props != null) {
                    geneFeatureTrack.setProperties(parser.getTrackProperties());
                }
                geneFeatureTrack.setUrl(annotationURL);
            }
        }
        return geneFeatureTrack;
    }

    /**
     * Create an annotation track for the genome from a supplied list of features
     *
     * @param genome
     * @param features
     */
    public FeatureTrack createGeneTrack(Genome genome, List<org.broad.tribble.Feature> features) {

        FeatureDB.clearFeatures();
        FeatureTrack geneFeatureTrack = null;
        String name = "Annotations";

        String id = genome.getId() + "_genes";
        geneFeatureTrack = new FeatureTrack(id, name, new FeatureCollectionSource(features, genome));
        geneFeatureTrack.setMinimumHeight(5);
        geneFeatureTrack.setHeight(35);
        //geneFeatureTrack.setRendererClass(GeneRenderer.class);
        geneFeatureTrack.setColor(Color.BLUE.darker());

        return geneFeatureTrack;
    }

    public void excludedUrl(String location) {
        excludedArchivesUrls.add(location);
    }

}
