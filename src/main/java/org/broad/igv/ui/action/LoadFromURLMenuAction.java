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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeDownloadUtils;
import org.broad.igv.feature.genome.load.HubGenomeLoader;
import org.broad.igv.logging.*;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.session.SessionReader;
import org.broad.igv.ucsc.hub.Hub;
import org.broad.igv.ucsc.hub.HubParser;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMenuBar;
import org.broad.igv.ui.util.LoadFromURLDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class LoadFromURLMenuAction extends MenuAction {

    public static final String LOAD_TRACKS_FROM_URL = "Load Tracks from URL...";
    public static final String LOAD_SESSION_FROM_URL = "Load Session from URL...";
    public static final String LOAD_GENOME_FROM_URL = "Load Genome from URL...";
    public static final String LOAD_FROM_HTSGET = "Load from htsget Server...";
    public static final String LOAD_TRACKHUB = "Load Track Hub...";

    private static final Logger log = LogManager.getLogger(LoadFromURLMenuAction.class);
    private final IGV igv;

    public LoadFromURLMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        JPanel ta = new JPanel();
        ta.setPreferredSize(new Dimension(600, 20));
        String command = e.getActionCommand();
        boolean isHtsGet = command.equalsIgnoreCase(LOAD_FROM_HTSGET);
        if (command.equalsIgnoreCase(LOAD_TRACKS_FROM_URL) || isHtsGet) {

            LoadFromURLDialog dlg = new LoadFromURLDialog(IGV.getInstance().getMainFrame(), isHtsGet);
            dlg.setVisible(true);

            if (!dlg.isCanceled()) {
                loadUrls(dlg.getFileURLs(), dlg.getIndexURLs(), isHtsGet);
            }
        } else if ((command.equalsIgnoreCase(LOAD_SESSION_FROM_URL))) {

            final String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta, "Enter URL to .xml session file", JOptionPane.QUESTION_MESSAGE);

            if (url != null && !url.trim().isBlank()) {
                try {
                    LongRunningTask.submit(() -> this.igv.loadSession(url.trim(), null));
                } catch (Exception ex) {
                    MessageUtils.showMessage("Error loading session: " + ex.getMessage());
                }
            }
        } else if ((command.equalsIgnoreCase(LOAD_GENOME_FROM_URL))) {

            String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta, "Enter URL to .json, hub.txt, or FASTA file", JOptionPane.QUESTION_MESSAGE);

            loadGenomeFromUrl(url);

        } else if ((command.equalsIgnoreCase(LOAD_TRACKHUB))) {
            loadTrackHub(ta);
        }
    }

    private void loadUrls(List<String> inputs, List<String> indexes, boolean isHtsGet) {

        if (inputs.size() == 1 && HubGenomeLoader.isHubURL(inputs.getFirst())) {

            LongRunningTask.submit(() -> {
                try {
                    Genome genome = GenomeManager.getInstance().getCurrentGenome();
                    String id = genome != null ? genome.getUCSCId() : null;

                    Hub hub = HubParser.loadHub(inputs.getFirst(), id);
                    if (hub.isAssemblyHub() && (genome == null || !hub.getGenomeConfig().getUcscID().equals(id))) {
                        HubGenomeLoader.loadAssemblyHub(hub);
                    } else if (genome != null) {
                        SelectHubTracksAction.selectAndLoadTracks(hub);
                        genome.addTrackHub(hub);
                        IGVMenuBar.getInstance().updateTracksMenu(genome);
                        GenomeDownloadUtils.saveLocalGenome(genome.getConfig());
                    }

                } catch (IOException ex) {
                    log.error("Error loading tack hub", ex);
                    MessageUtils.showMessage("Error loading track hub: " + ex.getMessage());

                }
            });
        } else if (inputs.size() == 1 && SessionReader.isSessionFile(inputs.getFirst())) {
            // Session URL
            String url = inputs.getFirst();

            try {
                LongRunningTask.submit(() -> this.igv.loadSession(url, null));
            } catch (Exception ex) {
                MessageUtils.showMessage("Error loading url: " + url + " (" + ex + ")");
            }
        } else {
            if (!indexes.isEmpty() && indexes.size() != inputs.size()) {
                throw new RuntimeException("The number of Index URLs must equal the number of File URLs");
            }

            List<ResourceLocator> locators = getResourceLocators(inputs, indexes, isHtsGet);
            igv.addToRecentUrls(locators);
            igv.loadTracks(locators);
        }
    }

    // Note: this is not currently used
    private static void loadTrackHub(JPanel ta) {
        String urlOrAccension = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta, "Enter GCA or GCF accession, or URL to hub.txt file", JOptionPane.QUESTION_MESSAGE);

        if (urlOrAccension == null) return;
        urlOrAccension = urlOrAccension.trim();
        final String url;
        if (urlOrAccension.startsWith("GC")) {
            url = HubGenomeLoader.convertToHubURL(urlOrAccension);
            if (!FileUtils.resourceExists(url)) {
                MessageUtils.showMessage("Unrecognized hub identifier: " + urlOrAccension);
            }
        } else {
            url = urlOrAccension;
        }

        loadGenomeFromUrl(url);
    }

    private static void loadGenomeFromUrl(String url) {
        if (url != null && !url.isBlank()) {
            url = url.trim();
            try {
                if (isHubURL(url)) {
                    HubGenomeLoader.loadGenome(url);
                } else {
                    GenomeManager.getInstance().loadGenome(url);
                }
            } catch (Exception e) {
                MessageUtils.showMessage("Error loading genome: " + e.getMessage());
            }
        }
    }

    private static List<ResourceLocator> getResourceLocators(List<String> inputs, List<String> indexes, boolean isHtsGet) {
        List<ResourceLocator> locators = new ArrayList<>();
        for (int i = 0; i < inputs.size(); i++) {
            final String url = inputs.get(i);
            final ResourceLocator rl = new ResourceLocator(url.trim());
            if (!indexes.isEmpty()) {
                final String indexUrl = indexes.get(i);
                rl.setIndexPath(indexUrl);
            }
            rl.setHtsget(isHtsGet);
            locators.add(rl);
        }
        return locators;
    }

    /**
     * Somewhat crude test for a hub url
     *
     * @param input
     * @return
     */
    private static boolean isHubURL(String input) {
        return input.endsWith("/hub.txt");
    }

}

