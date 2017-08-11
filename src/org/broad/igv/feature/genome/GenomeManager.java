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


import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.event.GenomeChangeEvent;
import org.broad.igv.event.GenomeResetEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.fasta.FastaBlockCompressedSequence;
import org.broad.igv.feature.genome.fasta.FastaDirectorySequence;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.feature.genome.fasta.FastaUtils;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.commandbar.GenomeListManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.*;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.ui.util.download.Downloader;
import org.broad.igv.util.*;

import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.SocketException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.*;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;

/**
 * @author jrobinso
 */
public class GenomeManager {

    final static String GENOME_ARCHIVE_VERSION_KEY = "version";
    final static String GENOME_ARCHIVE_PROPERTY_FILE_NAME = "property.txt";
    final static String GENOME_ARCHIVE_ID_KEY = "id";
    final static String GENOME_ARCHIVE_NAME_KEY = "name";
    final static String GENOME_ORDERED_KEY = "ordered";
    final static String GENOME_GENETRACK_NAME = "geneTrackName";
    final static String GENOME_URL_KEY = "url";
    final static String GENOME_ARCHIVE_CYTOBAND_FILE_KEY = "cytobandFile";
    final static String GENOME_ARCHIVE_GENE_FILE_KEY = "geneFile";
    final static String GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY = "sequenceLocation";
    final static String COMPRESSED_SEQUENCE_PATH = "compressedSequencePath";

    /**
     * Whether the sequenceLocation has been modified from the version of the .genome
     * file on the server
     */
    public static final String GENOME_ARCHIVE_CUSTOM_SEQUENCE_LOCATION_KEY = "customSequenceLocation";
    public static final String GENOME_CHR_ALIAS_FILE_KEY = "chrAliasFile";
    public static final String SEQUENCE_MAP_FILE = "sequenceMap.txt";
    private static Logger log = Logger.getLogger(GenomeManager.class);

    private static final String ACT_USER_DEFINED_GENOME_LIST_FILE = "user-defined-genomes.txt";
    public static final String TEST_USER_DEFINED_GENOME_LIST_FILE = "test-user-defined-genomes.txt";

    public static final GenomeListItem DEFAULT_GENOME = new GenomeListItem("Human hg19", "http://s3.amazonaws.com/igv.broadinstitute.org/genomes/hg19.genome", "hg19");

    private GenomeListManager genomeListManager;

//    public static String getUserDefinedGenomeListFile() {
//        if (Globals.isTesting()) {
//            return TEST_USER_DEFINED_GENOME_LIST_FILE;
//        } else {
//            return ACT_USER_DEFINED_GENOME_LIST_FILE;
//        }
//
//    }

    private static GenomeManager theInstance;

    private Genome currentGenome;

    boolean serverGenomeListUnreachable = false;

    private Map<String, File> localSequenceMap;

    /**
     * Map from genomeID -> GenomeListItem
     * ID comparison will be case insensitive
     */

    public synchronized static GenomeManager getInstance() {
        if (theInstance == null) {
            theInstance = new GenomeManager();
        }
        return theInstance;
    }

    private GenomeManager() {
        genomeListManager = GenomeListManager.getInstance();
        localSequenceMap = loadSequenceMap();
    }

    public void setCurrentGenome(Genome genome) {
        if (genome != null) {
            PreferencesManager.getPreferences().setLastGenome(genome.getId());
        }
        this.currentGenome = genome;
        if (genome != null) {
            if (IGV.hasInstance()) {
                IGV.getInstance().getSession().clearHistory();
                FrameManager.getDefaultFrame().setChromosomeName(genome.getHomeChromosome(), true);
            }
            IGVEventBus.getInstance().post(new GenomeChangeEvent(genome));
        }
    }

    public boolean isServerGenomeListUnreachable() {
        return serverGenomeListUnreachable;
    }

    public void loadGenomeById(String genomeId) throws IOException {

        final Genome currentGenome = getCurrentGenome();
        if (currentGenome != null && genomeId.equals(currentGenome.getId())) {
            return; // Already loaded
        }

        if (org.broad.igv.util.ParsingUtils.pathExists(genomeId)) {
            loadGenome(genomeId, null);

        } else {

            final ProgressMonitor[] monitor = {new ProgressMonitor()};
            final ProgressBar.ProgressDialog[] progressDialog = new ProgressBar.ProgressDialog[1];
            UIUtilities.invokeAndWaitOnEventThread(() -> {
                progressDialog[0] = ProgressBar.showProgressDialog(IGV.getMainFrame(), "Loading Genome...", monitor[0], false);
            });

            try {
                GenomeListItem item = genomeListManager.getGenomeListItem(genomeId);
                loadGenome(item.getPath(), monitor[0]);
            } finally {
                UIUtilities.invokeOnEventThread(() -> {
                    progressDialog[0].setVisible(false);
                });
            }


        }
    }

    public Genome loadGenome(String genomePath, ProgressMonitor monitor) throws IOException {

        try {
            log.info("Loading genome: " + genomePath);

            Genome newGenome = null;

            if (monitor != null) {
                UIUtilities.invokeAndWaitOnEventThread(() -> monitor.fireProgress(25));
            }

            // Clear Feature DB
            FeatureDB.clearFeatures();

            String altGenomePath;
            if (genomePath.endsWith(".genome")) {
                File archiveFile = getArchiveFile(genomePath);

                if(!archiveFile.exists()) {
                    return null;    // Happens if genome download was canceled.
                }

                altGenomePath = archiveFile.getAbsolutePath();
                newGenome = loadDotGenomeFile(archiveFile);
            } else if (genomePath.endsWith(".gbk") || genomePath.endsWith(".gb")) {
                altGenomePath = genomePath;
                newGenome = loadGenbankFile(genomePath);
            } else if (genomePath.endsWith(".chrom.sizes")) {
                altGenomePath = genomePath;
                newGenome = loadChromSizes(genomePath);
            } else {
                // Assume a fasta file
                altGenomePath = genomePath;
                if (genomePath.endsWith(Globals.GZIP_FILE_EXTENSION)) {

                    String gziPath = genomePath + ".gzi";
                    String faiPath = genomePath + ".fai";
                    if (!(FileUtils.resourceExists(gziPath) && FileUtils.resourceExists(faiPath))) {
                        throw new GenomeException("IGV cannot readed gzipped fasta files.");

                    }
                }
                if (!FileUtils.isRemote(genomePath)) {
                    if (!(new File(genomePath)).exists()) {
                        throw new GenomeException("Cannot locate genome: " + genomePath);
                    }
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
                monitor.fireProgress(25);
            }

            if (IGV.hasInstance()) IGV.getInstance().resetSession(null);


            GenomeListItem genomeListItem = new GenomeListItem(newGenome.getDisplayName(), altGenomePath, newGenome.getId());
            final Set<String> serverGenomeIDs = genomeListManager.getServerGenomeIDs();

            boolean userDefined = !serverGenomeIDs.contains(newGenome.getId());
            genomeListManager.addGenomeItem(genomeListItem, userDefined);


            setCurrentGenome(newGenome);

            if (IGV.hasInstance()) {
                FeatureTrack geneFeatureTrack = newGenome.getGeneTrack();
                IGV.getInstance().setGenomeTracks(geneFeatureTrack);
            }


            log.info("Genome loaded.  id= " + newGenome.getId());

            return currentGenome;

        } catch (SocketException e) {
            throw new RuntimeException("Server connection error", e);
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


        GenomeListItem item = genomeListManager.buildItemFromPath(fastaPath);
        if (item == null) {
            throw new IOException(fastaPath + " does not exist, could not load genome");
        }

        FastaIndexedSequence fastaSequence = fastaPath.endsWith(".gz") ?
                new FastaBlockCompressedSequence(fastaPath) :
                new FastaIndexedSequence(fastaPath);
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
     * Create a genome from a ".genome" file.  In addition to the reference sequence .genome files can optionally
     * specify cytobands and annotations.
     */
    private Genome loadDotGenomeFile(File archiveFile) throws IOException {

        Genome newGenome;

        GenomeDescriptor genomeDescriptor = parseGenomeArchiveFile(archiveFile);

        final String id = genomeDescriptor.getId();
        final String displayName = genomeDescriptor.getName();

        boolean isFasta = genomeDescriptor.isFasta();
        String[] fastaFiles = genomeDescriptor.getFastaFileNames();

        LinkedHashMap<String, List<Cytoband>> cytobandMap = null;
        if (genomeDescriptor.hasCytobands()) {
            cytobandMap = loadCytobandFile(genomeDescriptor);
        }

        String sequencePath = localSequenceMap.containsKey(genomeDescriptor.getId()) ?
                loadSequenceMap().get(genomeDescriptor.getId()).getAbsolutePath() :
                genomeDescriptor.getSequencePath();

        // Convert legacy "local fasta" .genome file
        if (genomeDescriptor.hasCustomSequenceLocation()) {
            String localPath = genomeDescriptor.getSequencePath();
            addLocalFasta(genomeDescriptor.getId(), new File(localPath));
        }

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

            if (sequencePath.endsWith(".gz")) {
                FastaBlockCompressedSequence fastaSequence = new FastaBlockCompressedSequence(sequencePath);
                sequence = new SequenceWrapper(fastaSequence);
            } else {
                FastaIndexedSequence fastaSequence = new FastaIndexedSequence(sequencePath);
                sequence = new SequenceWrapper(fastaSequence);
            }
            chromosOrdered = true;
        }

        newGenome = new Genome(id, displayName, sequence, chromosOrdered, genomeDescriptor);

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
    public void refreshCache(File cachedFile, URL genomeArchiveURL) {
        // Look in cache first
        try {
            if (cachedFile.exists()) {

                //File sizes won't be the same if the local version has a different sequence location
                boolean remoteModfied = HttpUtils.getInstance().remoteIsNewer(cachedFile, genomeArchiveURL);

                // Force an update of cached genome if file length does not equal remote content length
                boolean forceUpdate = remoteModfied &&
                        PreferencesManager.getPreferences().getAsBoolean(Constants.AUTO_UPDATE_GENOMES);

                if (forceUpdate) {

                    log.info("Refreshing genome: " + genomeArchiveURL.toString());

                    Downloader.download(genomeArchiveURL, cachedFile, IGV.getMainFrame());

                }
            } else {
                // Copy file directly from the server to local cache.
                Frame parent = IGV.hasInstance() ? IGV.getMainFrame() : null;
                Downloader.download(genomeArchiveURL, cachedFile, parent);
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

                if (zipEntryName.equalsIgnoreCase(GENOME_ARCHIVE_PROPERTY_FILE_NAME)) {
                    InputStream inputStream = zipFile.getInputStream(zipEntry);
                    Properties properties = new Properties();
                    properties.load(inputStream);

                    String cytobandZipEntryName = properties.getProperty(GENOME_ARCHIVE_CYTOBAND_FILE_KEY);
                    String geneFileName = properties.getProperty(GENOME_ARCHIVE_GENE_FILE_KEY);
                    String chrAliasFileName = properties.getProperty(GENOME_CHR_ALIAS_FILE_KEY);
                    String sequencePath = properties.getProperty(GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY);
                    String compressedSequencePath = properties.getProperty(COMPRESSED_SEQUENCE_PATH);

                    if ((sequencePath != null) && !HttpUtils.isRemoteURL(sequencePath)) {
                        sequencePath = getFullPath(f, sequencePath);
                    }

                    if ((compressedSequencePath != null) && !HttpUtils.isRemoteURL(compressedSequencePath)) {
                        compressedSequencePath = getFullPath(f, sequencePath);
                    }

                    boolean chrNamesAltered = parseBooleanPropertySafe(properties, "filenamesAltered");
                    boolean fasta = parseBooleanPropertySafe(properties, "fasta");
                    boolean fastaDirectory = parseBooleanPropertySafe(properties, "fastaDirectory");
                    boolean chromosomesAreOrdered = parseBooleanPropertySafe(properties, GENOME_ORDERED_KEY);
                    boolean hasCustomSequenceLocation = parseBooleanPropertySafe(properties, GENOME_ARCHIVE_CUSTOM_SEQUENCE_LOCATION_KEY);


                    String fastaFileNameString = properties.getProperty("fastaFiles");
                    String url = properties.getProperty(GENOME_URL_KEY);


                    // The new descriptor
                    genomeDescriptor = new GenomeZipDescriptor(
                            properties.getProperty(GENOME_ARCHIVE_NAME_KEY),
                            chrNamesAltered,
                            properties.getProperty(GENOME_ARCHIVE_ID_KEY),
                            cytobandZipEntryName,
                            geneFileName,
                            chrAliasFileName,
                            properties.getProperty(GENOME_GENETRACK_NAME, "Gene"),
                            sequencePath,
                            hasCustomSequenceLocation,
                            compressedSequencePath,
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

    private static String getFullPath(File f, String sequencePath) throws IOException {
        File sequenceFolder = null;
        sequenceFolder = new File(sequencePath);
        boolean isAbsolutePath = sequenceFolder.isAbsolute() ||
                sequencePath.startsWith("/") || sequencePath.startsWith("\\");
        if (!isAbsolutePath) {
            sequenceFolder = new File(f.getParent(), sequencePath);
        }
        sequencePath = sequenceFolder.getCanonicalPath();
        sequencePath.replace('\\', '/');
        return sequencePath;
    }

    private static boolean parseBooleanPropertySafe(Properties properties, String key) {
        String propertyString = properties.getProperty(key);
        return Boolean.parseBoolean(propertyString);
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
                                       javax.swing.ProgressMonitor monitor)
            throws IOException {

        File refFlatFile = null;
        File cytobandFile = null;
        File chrAliasFile = null;

        if (genomeFile != null) {
            PreferencesManager.getPreferences().setLastGenomeImportDirectory(genomeFile.getParentFile());
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

        if (monitor != null) monitor.setProgress(25);

        (new GenomeImporter()).createGenomeArchive(genomeFile, genomeId,
                genomeDisplayName, fastaFileName, refFlatFile, cytobandFile, chrAliasFile);

        if (monitor != null) monitor.setProgress(75);

        GenomeListItem newItem = new GenomeListItem(genomeDisplayName, genomeFile.getAbsolutePath(), genomeId);
        genomeListManager.addGenomeItem(newItem, true);

        if (monitor != null) monitor.setProgress(100);

        return newItem;

    }

    /**
     * Specific to Broad Amazon servers -- use S3 downwload rather than cloudfront
     *
     * @param path
     * @return
     */
    private String convertToS3(String path) {
        if (path.startsWith("http://igvdata") || path.startsWith("https://igvdata")) {
            return path.replaceFirst("igvdata", "igv");
        } else {
            return path;
        }
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


    public boolean downloadGenomes(final List<GenomeListItem> addValuesList, boolean downloadSequence) {

        boolean success = false;

        for (GenomeListItem item : addValuesList) {
            try {


                File archiveFile = getArchiveFile(item.getPath());                  // Has side affect of downloading .genome file

                if (downloadSequence && item.getPath().endsWith(".genome")) {

                    GenomeDescriptor genomeDescriptor = parseGenomeArchiveFile(archiveFile);

                    if (genomeDescriptor.isFasta()) {
                        String fastaPath = genomeDescriptor.getSequencePath();
                        File localFile = downloadFasta(fastaPath);
                        if(localFile != null) {
                            success = true;
                            addLocalFasta(item.getId(), localFile);
                        }

                    } else {
                        MessageUtils.showMessage("Could not download sequence for: " + genomeDescriptor.getName());
                    }
                }

                if(success) {
                    genomeListManager.addGenomeItem(item, false);
                }
            } catch (Exception e) {
                log.error("Fasta file unavailable for " + item.getDisplayableName());
            }
        }

        if(success) {
            IGVEventBus.getInstance().post(new GenomeResetEvent());
        }

        return success;

    }


    /**
     * Download a fasta file and associated index files.
     *
     * @throws IOException
     */
    File downloadFasta(String fastaPath) throws IOException {

        File defaultDir = DirectoryManager.getFastaCacheDirectory();
        File targetDir = defaultDir;
        //File targetDir = FileDialogUtils.chooseDirectory("Select directory for sequence", defaultDir);
        if (targetDir == null) {
            targetDir = defaultDir;
        }

        String filename = Utilities.getFileNameFromURL(fastaPath);

        File localFile = new File(targetDir, filename);
        boolean downloaded = Downloader.download(new URL(fastaPath), localFile, IGV.getMainFrame());

        if(downloaded) {
            URL indexUrl = new URL(fastaPath + ".fai");
            File localIndexFile = new File(targetDir, filename + ".fai");
            downloaded = Downloader.download(indexUrl, localIndexFile, IGV.getMainFrame());
        }

        if(downloaded) {

            if (fastaPath.endsWith(".gz")) {
                URL gziUrl = new URL(fastaPath + ".gzi");
                File localGziPath = new File(targetDir, filename + ".gzi");
                downloaded = Downloader.download(gziUrl, localGziPath, IGV.getMainFrame());
            }
        }

        return downloaded ? localFile : null;
    }


    private Map<String, File> loadSequenceMap() {

        File sequenceFile = new File(DirectoryManager.getGenomeCacheDirectory(), SEQUENCE_MAP_FILE);

        localSequenceMap = new HashMap<>();

        if (sequenceFile.exists()) {
            BufferedReader br = null;

            try {
                br = new BufferedReader(new FileReader(sequenceFile));
                String nextLine;
                while ((nextLine = br.readLine()) != null) {
                    String[] tokens = nextLine.split("\t");
                    if (tokens.length > 1) {
                        File file = new File(tokens[1]);
                        if (file.exists()) {
                            localSequenceMap.put(tokens[0], file);
                        } else {
                            log.info("Sequence file not found: " + file.getAbsolutePath());
                        }
                    }
                }
            } catch (IOException e) {
                log.error("Error loading sequence map file", e);
            } finally {
                if (br != null) try {
                    br.close();
                } catch (IOException e) {
                    log.error("Error closing sequenceMap file", e);
                }
            }
        }

        return localSequenceMap;
    }


    public File getLocalFasta(String id) {
        return localSequenceMap.get(id);
    }

    public void removeLocalFasta(String id) {
        localSequenceMap.remove(id);
        updateSequenceMapFile();
    }

    private void addLocalFasta(String id, File localFile) {
        localSequenceMap.put(id, localFile);
        updateSequenceMapFile();

    }


    private void updateSequenceMapFile() {

        PrintWriter pw = null;

        try {
            File sequenceFile = new File(DirectoryManager.getGenomeCacheDirectory(), SEQUENCE_MAP_FILE);
            pw = new PrintWriter(new BufferedWriter(new FileWriter(sequenceFile)));

            for (Map.Entry<String, File> entry : localSequenceMap.entrySet()) {
                pw.println(entry.getKey() + "\t" + entry.getValue());
            }
        } catch (IOException e) {
            log.error("Error writing sequence map", e);
        } finally {
            if (pw != null) pw.close();
        }
    }

}
