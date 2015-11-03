/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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

import com.google.common.collect.Iterables;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.*;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMenuBar;
import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.*;
import org.broad.igv.util.collections.CI;
import htsjdk.tribble.util.ParsingUtils;

import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.SocketException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.*;
import java.util.List;
import java.util.zip.*;

/**
 * @author jrobinso
 */
public class GenomeManager {

    private static Logger log = Logger.getLogger(GenomeManager.class);

    private static final String ACT_USER_DEFINED_GENOME_LIST_FILE = "user-defined-genomes.txt";
    public static final String TEST_USER_DEFINED_GENOME_LIST_FILE = "test-user-defined-genomes.txt";

    public static String getUserDefinedGenomeListFile() {
        if (Globals.isTesting()) {
            return TEST_USER_DEFINED_GENOME_LIST_FILE;
        } else {
            return ACT_USER_DEFINED_GENOME_LIST_FILE;
        }

    }

    private static GenomeManager theInstance;

    private Genome currentGenome;

    private LinkedHashSet<GenomeListItem> userDefinedGenomeArchiveList;
    private List<GenomeListItem> serverGenomeArchiveList;
    private List<GenomeListItem> cachedGenomeArchiveList;
    private Set<String> excludedArchivesUrls = new HashSet();

    /**
     * Map from genomeID -> GenomeListItem
     * ID comparison will be case insensitive
     */
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

    public boolean isServerGenomeListUnreachable() {
        return serverGenomeListUnreachable;
    }

    public Genome loadGenome(String genomePath, ProgressMonitor monitor) throws IOException {
        return loadGenome(genomePath, monitor, true);
    }

    /**
     * Load a genome from the given path.  Could be a .genome, .gbk, chrom.sizes, or fasta file
     *
     * @param genomePath     File, http, or ftp path to the .genome or indexed fasta file
     * @param monitor        ProgressMonitor  Monitor object, can be null
     * @param addGenomeTrack Whether to add the genomeTrack to IGV
     * @return Genome
     * @throws FileNotFoundException
     */
    public Genome loadGenome(String genomePath, ProgressMonitor monitor, boolean addGenomeTrack)
            throws IOException {

        try {
            log.info("Loading genome: " + genomePath);

            Genome newGenome = null;

            if (monitor != null) {
                monitor.fireProgressChange(25);
            }

            // Clear Feature DB
            FeatureDB.clearFeatures();

            if (genomePath.endsWith(".genome")) {
                newGenome = loadDotGenomeFile(genomePath);
            } else if (genomePath.endsWith(".gbk") || genomePath.endsWith(".gb")) {
                newGenome = loadGenbankFile(genomePath);
            } else if (genomePath.endsWith(".chrom.sizes")) {
                newGenome = loadChromSizes(genomePath);
            } else {
                // Assume a fasta file
                if (genomePath.endsWith(Globals.GZIP_FILE_EXTENSION)) {
                    throw new GenomeException("IGV cannot readed gzipped fasta files.");
                }
                newGenome = loadFastaFile(genomePath);
            }

            // Load alias files from genome source directory, if any
            String aliasPath = FileUtils.getParent(genomePath) + "/" + newGenome.getId() + "_alias.tab";
            Collection<Collection<String>> aliases = loadChrAliases(aliasPath);
            if (aliases != null) newGenome.addChrAliases(aliases);


            // Load user-defined chr aliases, if any.  This is done last so they have priority
            aliasPath = (new File(DirectoryManager.getGenomeCacheDirectory(), newGenome.getId() + "_alias.tab")).getAbsolutePath();
            aliases = loadChrAliases(aliasPath);
            if (aliases != null) newGenome.addChrAliases(aliases);

            if (monitor != null) {
                monitor.fireProgressChange(25);
            }

            setCurrentGenome(newGenome);

            if (IGV.hasInstance() && !Globals.isHeadless() && addGenomeTrack) {
                FeatureTrack geneFeatureTrack = newGenome.getGeneTrack();
                IGV.getInstance().setGenomeTracks(geneFeatureTrack);
            }

            log.info("Genome loaded.  id= " + newGenome.getId());

            return currentGenome;

        } catch (SocketException e) {
            throw new GenomeServerException("Server connection error", e);
        }

    }

    /**
     * Define a minimal genome from a chrom.sizes file.  It is assumed (required) that the file follow the
     * UCSC naming convention  =>  [id].chrom.sizes
     *
     * @param genomePath
     * @return
     * @throws IOException
     */
    private Genome loadChromSizes(String genomePath) throws IOException {

        int firstPeriodIdx = genomePath.indexOf('.');
        String genomeId = genomePath.substring(0, firstPeriodIdx);
        List<Chromosome> chromosomes = ChromSizesParser.parse(genomePath);
        Genome newGenome = new Genome(genomeId, chromosomes);

        // Search for chr aliases

        setCurrentGenome(newGenome);
        return newGenome;

    }

    private Genome loadGenbankFile(String genomePath) throws IOException {
        Genome newGenome;
        GenbankParser genbankParser = new GenbankParser(genomePath);
        genbankParser.readFeatures(true);

        String name = genbankParser.getLocusName();
        String chr = genbankParser.getChr();

        if (!name.equals(chr)) {
            name = name + " (" + chr + ")";
        }

        byte[] seq = genbankParser.getSequence();
        Sequence sequence = new InMemorySequence(chr, seq);
        newGenome = new Genome(chr, name, sequence, true);

        String[] aliases = genbankParser.getAliases();
        if (aliases != null) {
            List<String> aliasList = new ArrayList<String>();
            aliasList.add(chr);
            for (String a : aliases) {
                aliasList.add(a);
            }
            newGenome.addChrAliases(Arrays.<Collection<String>>asList(aliasList));
        }


        setCurrentGenome(newGenome);

        if (IGV.hasInstance() && !Globals.isHeadless()) {
            FeatureTrack geneFeatureTrack = createGeneTrack(newGenome, genbankParser.getFeatures());
            newGenome.setGeneTrack(geneFeatureTrack);
        }

        FeatureDB.addFeatures(genbankParser.getFeatures(), newGenome);

        return newGenome;
    }

    /**
     * Create a Genome from a single fasta file.
     *
     * @param genomePath
     * @return
     * @throws IOException
     */
    private Genome loadFastaFile(String genomePath) throws IOException {
        Genome newGenome;// Assume its a fasta
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
        newGenome = new Genome(item.getId(), item.getDisplayableName(), sequence, true);
        setCurrentGenome(newGenome);
        return newGenome;
    }

    private Collection<Collection<String>> loadChrAliases(String path) {

        // String id = genome.getId();
        // File aliasFile = new File(DirectoryManager.getGenomeCacheDirectory(), id + "_alias.tab");
        File aliasFile = new File(path);

        if (aliasFile.exists()) {

            BufferedReader br = null;

            try {
                br = new BufferedReader(new FileReader(aliasFile));
                return loadChrAliases(br);
            } catch (IOException e) {
                log.error("Error loading chr alias table", e);
                if (!Globals.isHeadless())
                    MessageUtils.showMessage("<html>Error loading chromosome alias table.  Aliases will not be avaliable<br>" +
                            e.toString());
            } finally {
                if (br != null) {
                    try {
                        br.close();
                    } catch (IOException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }
                }
            }
        }
        return null;
    }

    /**
     * @param path
     * @return GenomeListItem representing this path, or null
     * if it's a local path which doesn't exist (we don't check for existence of remote file)
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
    private Genome loadDotGenomeFile(String genomePath) throws IOException {
        Genome newGenome;
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
        boolean chromosOrdered = false;
        if (sequencePath == null) {
            sequence = null;
        } else if (!isFasta) {
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
            chromosOrdered = true;
        }

        newGenome = new Genome(id, displayName, sequence, chromosOrdered);
        if (cytobandMap != null) {
            newGenome.setCytobands(cytobandMap);
        }

        Collection<Collection<String>> aliases = loadChrAliases(genomeDescriptor);
        if (aliases != null) {
            newGenome.addChrAliases(aliases);
        }

        InputStream geneStream = null;
        String geneFileName = genomeDescriptor.getGeneFileName();
        if (geneFileName != null) {
            try {
                geneStream = genomeDescriptor.getGeneStream();
                if (geneFileName.endsWith(".gbk")) {
                    GenbankParser genbankParser = new GenbankParser();
                    genbankParser.readFeatures(geneStream, false);
                    FeatureTrack geneFeatureTrack = createGeneTrack(newGenome, genbankParser.getFeatures());
                    newGenome.setGeneTrack(geneFeatureTrack);
                } else {
                    BufferedReader reader = new BufferedReader(new InputStreamReader(geneStream));
                    FeatureTrack geneFeatureTrack = createGeneTrack(newGenome, reader,
                            geneFileName, genomeDescriptor.getGeneTrackName(),
                            genomeDescriptor.getUrl());

                    newGenome.setGeneTrack(geneFeatureTrack);
                }
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

    private static Collection<Collection<String>> loadChrAliases(BufferedReader br) throws IOException {
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
     * @return The chromosome alias table, or null if none is defined.
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
     * Refresh a locally cached genome if appropriate (newer one on server, user set preference to enable it)
     * If it doesn't have a local cache, just downloaded
     * If the cached version has a custom sequence location, that is copied over to the downloaded version
     *
     * @param cachedFile
     * @param genomeArchiveURL
     * @throws IOException
     */
    void refreshCache(File cachedFile, URL genomeArchiveURL) {
        // Look in cache first
        try {
            if (cachedFile.exists()) {
                //Check to see cached version has a custom sequence
                GenomeDescriptor cachedDescriptor = parseGenomeArchiveFile(cachedFile);
                //File sizes won't be the same if the local version has a different sequence location
                boolean remoteModfied = HttpUtils.getInstance().remoteIsNewer(cachedFile, genomeArchiveURL,
                        !cachedDescriptor.hasCustomSequenceLocation());

                // Force an update of cached genome if file length does not equal remote content length
                boolean forceUpdate = remoteModfied &&
                        PreferenceManager.getInstance().getAsBoolean(PreferenceManager.AUTO_UPDATE_GENOMES);
                if (forceUpdate) {
                    log.info("Refreshing genome: " + genomeArchiveURL.toString());
                    File tmpFile = new File(cachedFile.getAbsolutePath() + ".tmp");
                    if (HttpUtils.getInstance().downloadFile(genomeArchiveURL.toExternalForm(), tmpFile).isSuccess()) {

                        tmpFile.deleteOnExit();
                        boolean success = true;

                        if (cachedDescriptor.hasCustomSequenceLocation()) {
                            success = rewriteSequenceLocation(tmpFile, cachedDescriptor.getSequenceLocation());
                        }

                        if (success) {
                            FileUtils.copyFile(tmpFile, cachedFile);
                        } else {
                            log.warn("Updating genome failed: " + genomeArchiveURL.toString());
                        }
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
    public static GenomeDescriptor parseGenomeArchiveFile(File f)
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
                    boolean hasCustomSequenceLocation = parseBooleanPropertySafe(properties, Globals.GENOME_ARCHIVE_CUSTOM_SEQUENCE_LOCATION_KEY);


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
                            hasCustomSequenceLocation,
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

    private static boolean parseBooleanPropertySafe(Properties properties, String key) {
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
            GenomeManager.getInstance().addGenomeItems(Arrays.asList(matchingItem), false);
        }
        return matchingItem != null;
    }


    /**
     * Calls {@link #getServerGenomeArchiveList(Set)} with default set of excluded URLs
     *
     * @return
     */
    public List<GenomeListItem> getServerGenomeArchiveList() {
        return getServerGenomeArchiveList(excludedArchivesUrls);
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

        if (IGVMenuBar.getInstance() != null) {
            IGVMenuBar.getInstance().notifyGenomeServerReachable(!serverGenomeListUnreachable);
        }
        return serverGenomeArchiveList;
    }

    /**
     * Searches through currently loaded GenomeListItems and returns
     * that with a matching ID. null if not found. To search through
     * all server and user defined genomes, use {@link #findGenomeListItemById(String)}
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

    static GenomeListItem searchGenomeList(String genomeId, Iterable<GenomeListItem> genomeList) {
        if (genomeList == null) return null;
        for (GenomeListItem item : genomeList) {
            if (item.getId().equals(genomeId)) {
                return item;
            }
        }
        return null;
    }

    public List<GenomeListItem> getGenomes() {
        return new ArrayList<GenomeListItem>(genomeItemMap.values());
    }

    /**
     * Completely rebuild the genome drop down info.
     * This is based on preferences only, does not contact server
     */
    public void buildGenomeItemList() {
        genomeItemMap.clear();
        Collection<GenomeListItem> tmpuserDefinedGenomeList = null;
        Collection<GenomeListItem> tmpcachedGenomeArchiveList = null;

        try {
            tmpcachedGenomeArchiveList = getCachedGenomeArchiveList();
            tmpuserDefinedGenomeList = getUserDefinedGenomeArchiveList();
        } catch (IOException e) {
            MessageUtils.showErrorMessage("Cannot access user defined genome archive list", e);
        }

        String[] genomeIdArray = PreferenceManager.getInstance().getGenomeIdDisplayList();

        if (genomeIdArray.length == 0) {
            genomeIdArray = new String[]{PreferenceManager.getInstance().getDefaultGenome(), "hg18"};
        }

        Iterable<GenomeListItem> combinedTmp = Iterables.concat(tmpuserDefinedGenomeList, tmpcachedGenomeArchiveList);
        addGenomesToMap(genomeIdArray, combinedTmp, genomeItemMap);

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
    private void addGenomesToMap(String[] keepGenomeIds, Iterable<GenomeListItem> genomeListItems, Map<String, GenomeListItem> genomeMap) {
        for (String id : keepGenomeIds) {
            GenomeListItem genomeListItem = searchGenomeList(id, genomeListItems);

            //If we didn't find the id, it may be a path
            if (genomeListItem == null) {
                genomeListItem = buildFromPath(id);
            }

            if (genomeListItem == null) {
                /**
                 * The {@code id} isn't in {@code genomeListItems},
                 * and it's not a remote path or an existing local path
                 * We assume it's the ID for a genome stored on the server
                 */
                genomeListItem = new GenomeListItem(id, null, id);
            }

            if (!genomeMap.containsKey(id)) {
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

            File listFile = new File(DirectoryManager.getGenomeCacheDirectory(), getUserDefinedGenomeListFile());

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
                        userDefinedGenomeArchiveList.add(item);
                    } catch (Exception e) {
                        log.error("Error updating user genome list line '" + nextLine + "'", e);
                    }
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
     * Delete .genome files from the cache directory
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
        addGenomeItem(newItem, true);
        return newItem;

    }


    /**
     * Downloads .genome file, INCLUDING fasta file it points to,
     * and rewrites the property.txt file to point to the local version
     *
     * @param srcPath
     * @param targetDir
     * @throws IOException
     */
    public RunnableResult downloadWholeGenome(String srcPath, File targetDir, Frame dialogsParent) throws IOException {

        //Whether the srcGenome is simply a fasta file
        boolean isFastaFile = FastaUtils.isFastaPath(srcPath);
        boolean showProgressDialog = !Globals.isHeadless() && !Globals.isBatch();
        String srcFileName = Utilities.getFileNameFromURL(srcPath);

        if (isFastaFile) {
            //Simple case, just need to copy fasta and create index
            return downloadFasta(srcPath, targetDir, srcFileName, dialogsParent);

        } else if (srcFileName.endsWith(Globals.GENOME_FILE_EXTENSION)) {
            //Most useful case of a .genome file, which contains nearly everything in it, except the
            //sequence which is a fasta file stored elsewhere
            String genomeName = Utilities.getFileNameFromURL(srcPath);
            File srcGenomeArchive = new File(targetDir, genomeName);
            RunnableResult genomeResult = HttpUtils.getInstance().downloadFile(srcPath, srcGenomeArchive);
            if (!genomeResult.isSuccess()) return genomeResult;

            GenomeDescriptor descriptor = parseGenomeArchiveFile(srcGenomeArchive);

            File targetGenomeFile = new File(targetDir, srcGenomeArchive.getName());
            String sequencePath = descriptor.getSequenceLocation();

            boolean seqIsFasta = FastaUtils.isFastaPath(sequencePath);
            if (!seqIsFasta) {
                String msg = ("This genome sequence is not available for download. \n" +
                        "Please contact the igv team at https://groups.google.com/forum/#!forum/igv-help for further assistance");
                MessageUtils.showMessage(msg);
                return RunnableResult.CANCELLED;
            }

            String sequenceFileName = Utilities.getFileNameFromURL(sequencePath);
            File localSequenceFile = new File(targetDir, sequenceFileName);

            // Copy file directly from the server to local area
            // Shows cancellable dialog
            RunnableResult fastaResult = downloadFasta(sequencePath, targetDir, sequenceFileName, dialogsParent);
            if (fastaResult.isSuccess()) {
                //Rewrite properties file to point to local fasta
                rewriteSequenceLocation(targetGenomeFile, localSequenceFile.getAbsolutePath());
            }
            return fastaResult;
        } else {
            throw new IllegalArgumentException("Unknown file type. Cannot download " + srcPath);
        }
    }


    /**
     * Download a fasta file. Also attempts to download the index, and create an index
     * if that download fails
     *
     * @param path     Full source path of fasta
     * @param targetDir     Destination directory for fasta file
     * @param targetName    Destination file name for fasta file
     * @param dialogsParent Parent of progress dialog shown. None shown if null
     * @return
     * @throws IOException
     */
    private RunnableResult downloadFasta(String path, File targetDir, String targetName, Frame dialogsParent) throws IOException {
        //Simple case, just need to copy fasta and create index
        File destFile = new File(targetDir, targetName);

        String fastaPath = convertToS3(path);

        //TODO PROMPT TO OVERWRITE IF FILE EXISTS
        URLDownloader urlDownloader = HttpUtils.getInstance().downloadFile(fastaPath, destFile, dialogsParent, "Downloading genome sequence");
        RunnableResult fastaResult = urlDownloader.getResult();
        //If not successful for whatever reason, we don't get the index
        if (!fastaResult.isSuccess()) return fastaResult;

        File destIndexFile = new File(destFile.getAbsolutePath() + ".fai");
        String srcIndexPath = fastaPath + ".fai";
        RunnableResult idxResult = null;
        if (ParsingUtils.resourceExists(srcIndexPath)) {
            idxResult = HttpUtils.getInstance().downloadFile(srcIndexPath, destIndexFile);
        }
        //If remote fasta doesn't exist or download failed, we create our own index
        if (idxResult == null || idxResult == RunnableResult.FAILURE) {
            FastaUtils.createIndexFile(destFile.getAbsolutePath(), destIndexFile.getAbsolutePath());
        }

        return fastaResult;
    }

    /**
     * Specific to Broad Amazon servers -- use S3 downwload rather than cloudfront
     * @param path
     * @return
     */
    private String convertToS3(String path) {
        if(path.startsWith("http://igvdata") || path.startsWith("https://igvdata")) {
            return path.replaceFirst("igvdata", "igv");
        }
        else {
            return path;
        }
    }

    /**
     * Rewrite the {@link Globals#GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY} property to equal
     * the specified {@code newSequencePath}. Works by creating a temp file and renaming
     *
     * @param targetFile      A .genome file, in zip format
     * @param newSequencePath
     * @return boolean indicating success or failure.
     * @throws IOException
     */
    static boolean rewriteSequenceLocation(File targetFile, String newSequencePath) throws IOException {

        ZipFile targetZipFile = new ZipFile(targetFile);
        boolean success = false;

        File tmpZipFile = File.createTempFile("tmpGenome", ".zip");
        ZipEntry propEntry = targetZipFile.getEntry(Globals.GENOME_ARCHIVE_PROPERTY_FILE_NAME);

        InputStream propertyInputStream = null;
        ZipOutputStream zipOutputStream = null;
        Properties inputProperties = new Properties();


        try {
            propertyInputStream = targetZipFile.getInputStream(propEntry);
            BufferedReader reader = new BufferedReader(new InputStreamReader(propertyInputStream));

            //Copy over property.txt, only replacing a few properties
            inputProperties.load(reader);
            inputProperties.put(Globals.GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY, newSequencePath);
            inputProperties.put(Globals.GENOME_ARCHIVE_CUSTOM_SEQUENCE_LOCATION_KEY, Boolean.TRUE.toString());


            ByteArrayOutputStream propertyBytes = new ByteArrayOutputStream();
            PrintWriter propertyFileWriter = new PrintWriter(new OutputStreamWriter(propertyBytes));

            inputProperties.store(propertyFileWriter, null);

            propertyFileWriter.flush();
            byte[] newPropertyBytes = propertyBytes.toByteArray();

            Enumeration<? extends ZipEntry> entries = targetZipFile.entries();
            zipOutputStream = new ZipOutputStream(new FileOutputStream(tmpZipFile));
            while (entries.hasMoreElements()) {
                ZipEntry curEntry = entries.nextElement();
                ZipEntry writeEntry = null;

                if (curEntry.getName().equals(Globals.GENOME_ARCHIVE_PROPERTY_FILE_NAME)) {
                    writeEntry = new ZipEntry(Globals.GENOME_ARCHIVE_PROPERTY_FILE_NAME);
                    writeEntry.setSize(newPropertyBytes.length);
                    zipOutputStream.putNextEntry(writeEntry);
                    zipOutputStream.write(newPropertyBytes);
                    continue;
                } else {
                    //Because the compressed size can vary,
                    //we generate a new ZipEntry and copy some attributes
                    writeEntry = new ZipEntry(curEntry.getName());
                    writeEntry.setSize(curEntry.getSize());
                    writeEntry.setComment(curEntry.getComment());
                    writeEntry.setTime(curEntry.getTime());
                }

                zipOutputStream.putNextEntry(writeEntry);
                InputStream tmpIS = null;
                try {
                    tmpIS = targetZipFile.getInputStream(writeEntry);
                    int bytes = IOUtils.copy(tmpIS, zipOutputStream);
                    log.debug(bytes + " bytes written to " + targetFile);
                } finally {
                    if (tmpIS != null) tmpIS.close();
                }

            }
        } catch (Exception e) {
            tmpZipFile.delete();
            throw new RuntimeException(e.getMessage(), e);
        } finally {
            if (propertyInputStream != null) propertyInputStream.close();
            if (zipOutputStream != null) {
                zipOutputStream.flush();
                zipOutputStream.finish();
                zipOutputStream.close();
            }
            zipOutputStream = null;
            System.gc();
            success = true;
        }

        //This is a hack. I don't know why it's necessary,
        //but for some reason the output zip file seems to be corrupt
        //at least when called from GenomeManager.refreshArchive
        try {
            Thread.sleep(1500);
        } catch (InterruptedException e) {
            //
        }

        //Rename tmp file
        if (success) {
            targetFile.delete();
            FileUtils.copyFile(tmpZipFile, targetFile);
            success = targetFile.exists();
            tmpZipFile.delete();
        }
        return success;
    }

    public String getGenomeId() {
        return currentGenome == null ? null : currentGenome.getId();
    }

    /**
     * IGV always has exactly 1 genome loaded at a time.
     * This returns the currently loaded genome
     *
     * @return
     * @api
     */
    public Genome getCurrentGenome() {
        return currentGenome;
    }

    public void addGenomeItems(Collection<GenomeListItem> genomeListItems, boolean userDefined) {
        for (GenomeListItem genomeListItem : genomeListItems) {
            genomeItemMap.put(genomeListItem.getId(), genomeListItem);
            if (userDefined) userDefinedGenomeArchiveList.add(genomeListItem);
        }
        PreferenceManager.getInstance().saveGenomeIdDisplayList(genomeItemMap.values());
        updateImportedGenomePropertyFile();
    }

    public void addGenomeItem(GenomeListItem genomeListItem, boolean userDefined) {
        genomeItemMap.put(genomeListItem.getId(), genomeListItem);
        PreferenceManager.getInstance().saveGenomeIdDisplayList(genomeItemMap.values());
        if (userDefined) userDefinedGenomeArchiveList.add(genomeListItem);
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
                parser = new GFFParser();
            } else {
                parser = AbstractFeatureParser.getInstanceFor(new ResourceLocator(geneFileName), genome);
            }
            if (parser == null) {
                MessageUtils.showMessage("ERROR: Unrecognized annotation file format: " + geneFileName +
                        "<br>Annotations for genome: " + genome.getId() + " will not be loaded.");
            } else {
                List<htsjdk.tribble.Feature> genes = parser.loadFeatures(reader, genome);
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
    public FeatureTrack createGeneTrack(Genome genome, List<htsjdk.tribble.Feature> features) {

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


    /**
     * Delete the specified .genome files and their sequences, only if they were downloaded from the
     * server. Doesn't touch user defined genomes
     *
     * @param removedValuesList
     */
    public void deleteDownloadedGenomes(List<GenomeListItem> removedValuesList) throws IOException {
        Collection<GenomeListItem> userDefinedGenomes = getUserDefinedGenomeArchiveList();
        for (GenomeListItem item : removedValuesList) {
            if (userDefinedGenomes.contains(item)) {
                continue;
            }

            String loc = item.getLocation();
            if (!HttpUtils.isRemoteURL(loc)) {
                File genFile = new File(loc);
                GenomeDescriptor descriptor = parseGenomeArchiveFile(genFile);
                if (!HttpUtils.isRemoteURL(descriptor.getSequenceLocation())) {
                    File seqFile = new File(descriptor.getSequenceLocation());
                    seqFile.delete();
                    File indexFile = new File(seqFile.getAbsolutePath() + ".fai");
                    indexFile.delete();
                }
                genFile.delete();
            }
        }
    }
}
