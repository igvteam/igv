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


import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.event.GenomeChangeEvent;
import org.broad.igv.event.GenomeResetEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.genome.load.ChromAliasParser;
import org.broad.igv.feature.genome.load.GenomeDescriptor;
import org.broad.igv.feature.genome.load.GenomeLoader;
import org.broad.igv.feature.genome.load.JsonGenomeLoader;
import org.broad.igv.jbrowse.CircularViewUtilities;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.PanelName;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.commandbar.GenomeListManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.*;
import org.broad.igv.ui.util.download.Downloader;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;

import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.SocketException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.*;
import java.util.List;

import static org.broad.igv.prefs.Constants.SHOW_SINGLE_TRACK_PANE_KEY;

/**
 * @author jrobinso
 */
public class GenomeManager {

    private static Logger log = LogManager.getLogger(GenomeManager.class);

    private static GenomeManager theInstance;

    private static GenomeListManager genomeListManager;

    private Genome currentGenome;


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
        GenomeLoader.localSequenceMap = GenomeLoader.loadSequenceMap();
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
    public static File getGenomeFile(String genomePath) throws MalformedURLException, UnsupportedEncodingException {
        File archiveFile;
        if (HttpUtils.isRemoteURL(genomePath.toLowerCase())) {
            // We need a local copy, as there is no http zip file reader
            URL genomeArchiveURL = HttpUtils.createURL(genomePath);
            final String tmp = URLDecoder.decode(genomeArchiveURL.getFile(), "UTF-8");
            String cachedFilename = Utilities.getFileNameFromURL(tmp);
            if (!DirectoryManager.getGenomeCacheDirectory().exists()) {
                DirectoryManager.getGenomeCacheDirectory().mkdir();
            }
            archiveFile = new File(DirectoryManager.getGenomeCacheDirectory(), cachedFilename);
            Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
            Downloader.download(genomeArchiveURL, archiveFile, parent);
        } else {
            archiveFile = new File(genomePath);
        }
        return archiveFile;
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
                IGVEventBus.getInstance().post(new GenomeChangeEvent(genome));
            }
        }
    }


    public void loadGenomeById(String genomeId) throws IOException {
        final Genome currentGenome = getCurrentGenome();
        if (currentGenome != null && genomeId.equals(currentGenome.getId())) {
            return; // Already loaded
        }

        String genomePath = null;
        if (org.broad.igv.util.ParsingUtils.fileExists(genomeId)) {
            genomePath = genomeId;
        } else {
            GenomeListItem item = genomeListManager.getGenomeListItem(genomeId);
            if (item == null) {
                MessageUtils.showMessage("Could not locate genome with ID: " + genomeId);
                return;
            } else {
                genomePath = item.getPath();
            }
        }

        loadGenome(genomePath); // monitor[0]);

    }


    /**
     * The main load method -- loads a genome from a file or url path.  Note this is a long running operation and
     * should not be done on the Swing event thread as it will block the UI.
     *
     * @param genomePath
     * @return
     * @throws IOException
     */
    public Genome loadGenome(String genomePath) throws IOException {

        WaitCursorManager.CursorToken cursorToken = null;
        try {
            log.info("Loading genome: " + genomePath);
            if (IGV.hasInstance()) {
                IGV.getInstance().setStatusBarMessage("<html><font color=blue>Loading genome</font></html>");
                cursorToken = WaitCursorManager.showWaitCursor();
            }

            // Clear Feature DB
            FeatureDB.clearFeatures();

            Genome newGenome = GenomeLoader.getLoader(genomePath).loadGenome();

            // Load user-defined chr aliases, if any.  This is done last so they have priority
            try {
                String aliasPath = (new File(DirectoryManager.getGenomeCacheDirectory(), newGenome.getId() + "_alias.tab")).getAbsolutePath();
                if (!(new File(aliasPath)).exists()) {
                    aliasPath = (new File(DirectoryManager.getGenomeCacheDirectory(), newGenome.getId() + "_alias.tab.txt")).getAbsolutePath();
                }
                if ((new File(aliasPath)).exists()) {
                    newGenome.addChrAliases(ChromAliasParser.loadChrAliases(aliasPath));
                }
            } catch (Exception e) {
                log.error("Failed to load user defined alias", e);
            }


            if (IGV.hasInstance()) {
                IGV.getInstance().resetSession(null);
            }

            GenomeListItem genomeListItem = new GenomeListItem(newGenome.getDisplayName(), genomePath, newGenome.getId());
            final Set<String> serverGenomeIDs = genomeListManager.getServerGenomeIDs();

            boolean userDefined = !serverGenomeIDs.contains(newGenome.getId());
            genomeListManager.addGenomeItem(genomeListItem, userDefined);

            setCurrentGenome(newGenome);

            // hasInstance() test needed for unit tests
            if (IGV.hasInstance()) {
                IGV.getInstance().goToLocus(newGenome.getHomeChromosome()); //  newGenome.getDefaultPos());
                loadGenomeAnnotations(newGenome);
                IGV.getInstance().resetFrames();
            }

            if (PreferencesManager.getPreferences().getAsBoolean(Constants.CIRC_VIEW_ENABLED) && CircularViewUtilities.ping()) {
                CircularViewUtilities.changeGenome(newGenome);
            }

            // log.warn("Genome loaded.  id= " + newGenome.getId());
            return currentGenome;

        } catch (SocketException e) {
            throw new RuntimeException("Server connection error", e);
        } finally {
            if (IGV.hasInstance()) {
                IGV.getInstance().setStatusBarMessage("");
                WaitCursorManager.removeWaitCursor(cursorToken);
            }
        }
    }

    /**
     * Load and initialize the track objects from the genome's track resource locators.  Does not add the tracks
     * to the IGV instance.
     *
     * @param genome
     */
    public void loadGenomeAnnotations(Genome genome) {
        restoreGenomeTracks(genome);
        IGV.getInstance().repaint();
    }

    /**
     * Add a genomes tracks to the IGV instance.
     *
     * @param genome
     */
    public void restoreGenomeTracks(Genome genome) {

        IGV.getInstance().setSequenceTrack();

        // Fetch the gene track, defined by .genome files.  In this format the genome data is encoded in the .genome file
        FeatureTrack geneFeatureTrack = genome.getGeneTrack();   // Can be null
        if (geneFeatureTrack != null) {
            PanelName panelName = PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY) ?
                    PanelName.DATA_PANEL : PanelName.FEATURE_PANEL;
            geneFeatureTrack.setAttributeValue(Globals.TRACK_NAME_ATTRIBUTE, geneFeatureTrack.getName());
            geneFeatureTrack.setAttributeValue(Globals.TRACK_DATA_FILE_ATTRIBUTE, "");
            geneFeatureTrack.setAttributeValue(Globals.TRACK_DATA_TYPE_ATTRIBUTE, geneFeatureTrack.getTrackType().toString());
            IGV.getInstance().addTracks(Arrays.asList(geneFeatureTrack), panelName);
        }

        List<ResourceLocator> resources = genome.getAnnotationResources();
        List<Track> annotationTracks = new ArrayList<>();
        if (resources != null) {
            for (ResourceLocator locator : resources) {
                try {
                    List<Track> tracks = IGV.getInstance().load(locator);
                    annotationTracks.addAll(tracks);
                } catch (DataLoadException e) {
                    log.error("Error loading genome annotations", e);
                }
            }
        }

        if (annotationTracks.size() > 0) {
            IGV.getInstance().addTracks(annotationTracks);
            for (Track track : annotationTracks) {
                ResourceLocator locator = track.getResourceLocator();
                String fn = "";
                if (locator != null) {
                    fn = locator.getPath();
                    int lastSlashIdx = fn.lastIndexOf("/");
                    if (lastSlashIdx < 0) {
                        lastSlashIdx = fn.lastIndexOf("\\");
                    }
                    if (lastSlashIdx > 0) {
                        fn = fn.substring(lastSlashIdx + 1);
                    }
                }
                track.setAttributeValue(Globals.TRACK_NAME_ATTRIBUTE, track.getName());
                track.setAttributeValue(Globals.TRACK_DATA_FILE_ATTRIBUTE, fn);
                track.setAttributeValue(Globals.TRACK_DATA_TYPE_ATTRIBUTE, track.getTrackType().toString());
            }
        }

        IGV.getInstance().revalidateTrackPanels();
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


    public boolean downloadGenome(GenomeListItem item, boolean downloadSequence) {

        boolean success;
        try {
            File genomeFile = getGenomeFile(item.getPath());                  // Has side affect of downloading .genome file
            if (downloadSequence) {
                String fastaPath = null;
                if (item.getPath().endsWith(".genome")) {
                    GenomeDescriptor genomeDescriptor = GenomeDescriptor.parseGenomeArchiveFile(genomeFile);
                    fastaPath = genomeDescriptor.getSequencePath();
                } else if (item.getPath().endsWith(".json")) {
                    JsonGenomeLoader.GenomeDescriptor desc = (new JsonGenomeLoader(item.getPath())).loadDescriptor();
                    fastaPath = desc.getFastaURL();
                }
                if (fastaPath != null && FileUtils.isRemote(fastaPath)) {
                    File localFile = downloadFasta(fastaPath);
                    if (localFile != null) {
                        addLocalFasta(item.getId(), localFile);
                    }
                }
            }

            success = true;

        } catch (Exception e) {
            success = false;
            MessageUtils.showErrorMessage("Error downloading genome", e);
            log.error("Error downloading genome " + item.getDisplayableName());
        }


        if (success) {
            genomeListManager.addGenomeItem(item, false);
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
        boolean downloaded = Downloader.download(HttpUtils.createURL(fastaPath), localFile, IGV.getInstance().getMainFrame());

        if (downloaded) {
            URL indexUrl = HttpUtils.createURL(fastaPath + ".fai");
            File localIndexFile = new File(targetDir, filename + ".fai");
            downloaded = Downloader.download(indexUrl, localIndexFile, IGV.getInstance().getMainFrame());
        }

        if (downloaded) {

            if (fastaPath.endsWith(".gz")) {
                URL gziUrl = HttpUtils.createURL(fastaPath + ".gzi");
                File localGziPath = new File(targetDir, filename + ".gzi");
                downloaded = Downloader.download(gziUrl, localGziPath, IGV.getInstance().getMainFrame());
            }
        }

        return downloaded ? localFile : null;
    }


    public static File getLocalFasta(String id) {
        return GenomeLoader.localSequenceMap.get(id);
    }

    public static void removeLocalFasta(String id) {
        GenomeLoader.localSequenceMap.remove(id);
        updateSequenceMapFile();
    }

    private static void addLocalFasta(String id, File localFile) {
        GenomeLoader.localSequenceMap.put(id, localFile);
        updateSequenceMapFile();
    }


    private static void updateSequenceMapFile() {

        PrintWriter pw = null;

        try {
            File sequenceFile = new File(DirectoryManager.getGenomeCacheDirectory(), GenomeDescriptor.SEQUENCE_MAP_FILE);
            pw = new PrintWriter(new BufferedWriter(new FileWriter(sequenceFile)));

            for (Map.Entry<String, File> entry : GenomeLoader.localSequenceMap.entrySet()) {
                pw.println(entry.getKey() + "\t" + entry.getValue());
            }
        } catch (IOException e) {
            log.error("Error writing sequence map", e);
        } finally {
            if (pw != null) pw.close();
        }
    }

}
