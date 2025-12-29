/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.action;

import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.feature.genome.load.HubGenomeLoader;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.session.SessionReader;
import org.igv.track.AttributeManager;
import org.igv.ucsc.hub.Hub;
import org.igv.ucsc.hub.HubParser;
import org.igv.ucsc.hub.HubRegistry;
import org.igv.ui.IGV;
import org.igv.ui.IGVMenuBar;
import org.igv.ui.util.LoadFromURLDialog;
import org.igv.ui.util.MessageUtils;
import org.igv.util.LongRunningTask;
import org.igv.util.ResourceLocator;

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
    public static final String LOAD_HUB_FROM_URL = "Add Track Hub from URL...";
    public static final String LOAD_SESSION_FROM_URL = "Load Session from URL...";
    public static final String LOAD_SAMPLEINFO_FROM_URL = "Load Sample Info from URL...";

    public static final String LOAD_GENOME_FROM_URL = "Load Genome from URL...";
    public static final String LOAD_FROM_HTSGET = "Load from htsget Server...";

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
                List<String> urls = dlg.getFileURLs();
                if (urls.size() == 1 && SessionReader.isSessionFile(urls.get(0))) {
                    LongRunningTask.submit(() -> this.igv.loadSession(urls.get(0), null));
                } else {
                    loadUrls(urls, dlg.getIndexURLs(), isHtsGet);
                }
            }
        } else if (command.equalsIgnoreCase(LOAD_SAMPLEINFO_FROM_URL)) {
            final String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta,
                    "Enter URL to sample info file", JOptionPane.PLAIN_MESSAGE);

            if (url != null && !url.trim().isBlank()) {
                try {
                    ResourceLocator locator = new ResourceLocator(url.trim());
                    LongRunningTask.submit(() -> {
                        AttributeManager.getInstance().loadSampleInfo(locator);
                        igv.revalidateTrackPanels();
                    });
                } catch (Exception ex) {
                    MessageUtils.showMessage("Error loading sample info: " + ex.getMessage());
                }
            }

        } else if ((command.equalsIgnoreCase(LOAD_SESSION_FROM_URL))) {

            final String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta,
                    "Enter URL to .xml session file", JOptionPane.PLAIN_MESSAGE);

            if (url != null && !url.trim().isBlank()) {
                try {
                    LongRunningTask.submit(() -> this.igv.loadSession(url.trim(), null));
                } catch (Exception ex) {
                    MessageUtils.showMessage("Error loading session: " + ex.getMessage());
                }
            }
        } else if ((command.equalsIgnoreCase(LOAD_GENOME_FROM_URL))) {

            String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta,
                    "Enter URL to .json, hub.txt, or FASTA file", JOptionPane.PLAIN_MESSAGE);
            if (url != null && !url.trim().isBlank()) {
                loadGenomeFromUrl(url.trim());
            }

        } else if ((command.equalsIgnoreCase(LOAD_HUB_FROM_URL))) {
            String url = JOptionPane.showInputDialog(IGV.getInstance().getMainFrame(), ta,
                    "Enter URL to a  hub.txt file", JOptionPane.PLAIN_MESSAGE);
            if (url != null && !url.trim().isBlank()) {
                loadTrackHub(url.trim());
            }
        }
    }

    private void loadUrls(List<String> inputs, List<String> indexes, boolean isHtsGet) {

        if (inputs.size() == 1 && HubGenomeLoader.isHubURL(inputs.getFirst())) {
            loadTrackHub(inputs.getFirst());

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

    /**
     * Load a track hub from a URL.  This method is called from the menu action and also from the
     * HubGenomeLoader when loading a genome.
     *
     * @param url
     */
    private static void loadTrackHub(final String url) {

        LongRunningTask.submit(() -> {
            try {
                Genome genome = GenomeManager.getInstance().getCurrentGenome();
                String id = genome != null ? genome.getUCSCId() : null;

                Hub hub = HubParser.loadHub(url);
                if (hub.isAssemblyHub() && (genome == null || !hub.getGenomeConfigs().get(0).getUcscID().equals(id))) {
                    HubGenomeLoader.loadAssemblyHub(hub);
                } else if (genome != null) {
                    //SelectHubTracksAction.selectAndLoadTracks(hub, id);
                    HubRegistry.addUserHub(hub);
                    IGVMenuBar.getInstance().updateMenus(genome);
                }

            } catch (IOException ex) {
                log.error("Error loading tack hub", ex);
                MessageUtils.showMessage("Error loading track hub: " + ex.getMessage());

            }
        });
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
            if (isHtsGet) {
                rl.setHtsget(true);  // Do not override setting if isHtsGet is false.  False means htsget status is unknown.
            }
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

