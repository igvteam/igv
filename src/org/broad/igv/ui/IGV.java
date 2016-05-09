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
 * IGV.java
 *
 * Represents an IGV instance.
 *
 * Note:  Currently, only one instance is allowed per JVM.
 *
 */
package org.broad.igv.ui;

import apple.dts.samplecode.osxadapter.OSXAdapter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.eventbus.Subscribe;
import com.jidesoft.swing.JideSplitPane;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.annotations.ForTesting;
import org.broad.igv.batch.BatchRunner;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.MaximumContigGenomeException;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.*;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.Preloader;
import org.broad.igv.peaks.PeakCommandBar;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionReader;
import org.broad.igv.session.UCSCSessionReader;
import org.broad.igv.session.IndexAwareSessionReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.dnd.GhostGlassPane;
import org.broad.igv.ui.event.*;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.*;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.*;
import org.broad.igv.variant.VariantTrack;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.io.*;
import java.lang.ref.SoftReference;
import java.lang.reflect.InvocationTargetException;
import java.net.NoRouteToHostException;
import java.net.URL;
import java.util.*;
import java.util.List;
import java.util.concurrent.Future;
import java.util.prefs.Preferences;

import static org.broad.igv.ui.WaitCursorManager.CursorToken;

/**
 * Represents an IGV instance, consisting of a main window and associated model.
 *
 * @author jrobinso
 */
public class IGV {

    private static Logger log = Logger.getLogger(IGV.class);
    private static IGV theInstance;

    // Window components
    private Frame mainFrame;
    private JRootPane rootPane;
    private IGVContentPane contentPane;
    private IGVMenuBar menuBar;

    private StatusWindow statusWindow;

    // Glass panes
    Component glassPane;
    GhostGlassPane dNdGlassPane;

    // Cursors
    public static Cursor fistCursor;
    public static Cursor zoomInCursor;
    public static Cursor zoomOutCursor;
    public static Cursor dragNDropCursor;

    //Session session;
    Session session;

    private GenomeManager genomeManager;

    /**
     * Attribute used to group tracks.  Normally "null".  Set from the "Tracks" menu.
     */
    private String groupByAttribute = null;


    private Map<String, List<Track>> overlayTracksMap = new HashMap();
    private Set<Track> overlaidTracks = new HashSet();

    public static final String DATA_PANEL_NAME = "DataPanel";
    public static final String FEATURE_PANEL_NAME = "FeaturePanel";


    // Misc state
    private LinkedList<String> recentSessionList = new LinkedList<String>();
    private boolean isExportingSnapshot = false;

    // Listeners
    Collection<SoftReference<TrackGroupEventListener>> groupListeners =
            Collections.synchronizedCollection(new ArrayList<SoftReference<TrackGroupEventListener>>());

    Collection<SoftReference<AlignmentTrackEventListener>> alignmentTrackListeners =
            Collections.synchronizedCollection(new ArrayList<SoftReference<AlignmentTrackEventListener>>());

    private List<JComponent> otherToolMenus = new ArrayList<JComponent>();
    private boolean rulerEnabled;

    /**
     * Add an entry to the "Tools" menu
     *
     * @param menu
     * @api
     */
    public void addOtherToolMenu(JComponent menu) {
        otherToolMenus.add(menu);
        if (menuBar != null) menuBar.refreshToolsMenu();
    }


    List<JComponent> getOtherToolMenus() {
        return otherToolMenus;
    }

    @Subscribe
    public void receiveViewChange(ViewChange.Result e) {
        repaintDataAndHeaderPanels();
        repaintStatusAndZoomSlider();
    }

    @Subscribe
    public void receiveViewChange(ViewChange.ChromosomeChangeResult e) {
        chromosomeChangeEvent(e.chrName, false);
    }

    public static IGV createInstance(Frame frame) {
        if (theInstance != null) {
            throw new RuntimeException("Only a single instance is allowed.");
        }
        theInstance = new IGV(frame);
        return theInstance;
    }

    public static IGV getInstance() {
        if (theInstance == null) {
            throw new RuntimeException("IGV has not been initialized.  Must call createInstance(Frame) first");
        }
        return theInstance;
    }

    @ForTesting
    static void destroyInstance() {
        IGVMenuBar.destroyInstance();
        theInstance = null;
    }

    public static boolean hasInstance() {
        return theInstance != null;
    }

    public static JRootPane getRootPane() {
        return getInstance().rootPane;
    }

    /**
     * The IGV GUI has one master frame containing all other elements.
     * This method returns that frame.
     *
     * @return
     * @api
     */
    public static Frame getMainFrame() {
        return getInstance().mainFrame;
    }


    /**
     * Creates new IGV
     */
    private IGV(Frame frame) {

        theInstance = this;

        genomeManager = GenomeManager.getInstance();

        mainFrame = frame;
        mainFrame.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent windowEvent) {
                windowCloseEvent();
            }

            @Override
            public void windowClosed(WindowEvent windowEvent) {
                windowCloseEvent();
            }

            private void windowCloseEvent() {
                PreferenceManager.getInstance().setApplicationFrameBounds(mainFrame.getBounds());
            }

            @Override
            public void windowLostFocus(WindowEvent windowEvent) {
                // Start & stop tooltip manager to force any tooltip windows to close.
                ToolTipManager.sharedInstance().setEnabled(false);
                ToolTipManager.sharedInstance().setEnabled(true);
                IGVPopupMenu.closeAll();
            }


            @Override
            public void windowDeactivated(WindowEvent windowEvent) {
                // Start & stop tooltip manager to force any tooltip windows to close.
                ToolTipManager.sharedInstance().setEnabled(false);
                ToolTipManager.sharedInstance().setEnabled(true);
                IGVPopupMenu.closeAll();
            }

            @Override
            public void windowActivated(WindowEvent windowEvent) {

            }

            @Override
            public void windowGainedFocus(WindowEvent windowEvent) {

            }
        });


        session = new Session(null);

        // Create cursors
        createHandCursor();
        createZoomCursors();
        createDragAndDropCursor();

        // Create components
        mainFrame.setTitle(UIConstants.APPLICATION_NAME);

        if (mainFrame instanceof JFrame) {
            JFrame jf = (JFrame) mainFrame;
            rootPane = jf.getRootPane();
        } else {
            rootPane = new JRootPane();
            mainFrame.add(rootPane);

        }
        contentPane = new IGVContentPane(this);
        menuBar = IGVMenuBar.createInstance(this);

        rootPane.setContentPane(contentPane);
        rootPane.setJMenuBar(menuBar);
        glassPane = rootPane.getGlassPane();
        glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        glassPane.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                super.mouseClicked(e);
                if (Globals.isDevelopment()) {
                    Component cButton = contentPane.getStatusBar().getCancelButton();
                    int tX = -(contentPane.getStatusBar().getX() + cButton.getX());
                    int tY = -(contentPane.getY() + contentPane.getStatusBar().getY());
                    e.translatePoint(tX, tY);
                    if (cButton.contains(e.getPoint())) {
                        contentPane.getStatusBar().getCancelButton().doClick();
                    }
                }
            }
        });
        dNdGlassPane = new GhostGlassPane();

        mainFrame.pack();

        // Set the application's previous location and size
        Dimension screenBounds = Toolkit.getDefaultToolkit().getScreenSize();
        Rectangle applicationBounds = PreferenceManager.getInstance().getApplicationFrameBounds();
        int state = PreferenceManager.getInstance().getAsInt(PreferenceManager.FRAME_STATE_KEY);

        if (applicationBounds == null || applicationBounds.getMaxX() > screenBounds.getWidth() ||
                applicationBounds.getMaxY() > screenBounds.getHeight()) {
            int width = Math.min(1150, (int) screenBounds.getWidth());
            int height = Math.min(800, (int) screenBounds.getHeight());
            applicationBounds = new Rectangle(0, 0, width, height);
        }

        //Certain components MUST be visible, so we set minimum size
        //{@link MainPanel#addDataPanel}
        mainFrame.setMinimumSize(new Dimension(300, 300));

        mainFrame.setExtendedState(state);
        mainFrame.setBounds(applicationBounds);
    }


    public void repaint() {
        mainFrame.repaint();
    }


    public GhostGlassPane getDnDGlassPane() {
        return dNdGlassPane;
    }

    public void startDnD() {
        rootPane.setGlassPane(dNdGlassPane);
        dNdGlassPane.setVisible(true);
    }

    public void endDnD() {
        rootPane.setGlassPane(glassPane);
        glassPane.setVisible(false);
    }

    public Dimension getPreferredSize() {
        return UIConstants.preferredSize;
    }


    public void addRegionOfInterest(RegionOfInterest roi) {
        session.addRegionOfInterestWithNoListeners(roi);
        RegionOfInterestPanel.setSelectedRegion(roi);
        doRefresh();
    }

    void beginROI(JButton button) {
        for (TrackPanel tp : getTrackPanels()) {
            TrackPanelScrollPane tsv = tp.getScrollPane();
            DataPanelContainer dpc = tsv.getDataPanel();
            for (Component c : dpc.getComponents()) {
                if (c instanceof DataPanel) {
                    DataPanel dp = (DataPanel) c;
                    RegionOfInterestTool regionOfInterestTool = new RegionOfInterestTool(dp, button);
                    dp.setCurrentTool(regionOfInterestTool);
                }
            }
        }


    }

    public void endROI() {
        for (TrackPanel tp : getTrackPanels()) {
            DataPanelContainer dp = tp.getScrollPane().getDataPanel();
            dp.setCurrentTool(null);
        }

    }

    private void chromosomeChangeEvent(String chrName, boolean updateCommandBar) {
        contentPane.chromosomeChanged(chrName);
        repaintDataAndHeaderPanels(updateCommandBar);
        contentPane.getCommandBar().updateComponentStates();
    }

    public void repaintStatusAndZoomSlider() {
        contentPane.getCommandBar().updateComponentStates();
        contentPane.getCommandBar().repaint();
    }

    /**
     * Repaints dataAndHeaderPanels as well as
     * zoom controls IFF IGV has instance && not headless.
     * Mostly use for testing
     */
    public static void repaintPanelsHeadlessSafe() {
        if (IGV.hasInstance() && !Globals.isHeadless()) {
            IGV.getInstance().repaintDataAndHeaderPanels();
            IGV.getInstance().repaintStatusAndZoomSlider();
        }
    }

    /**
     * Repaint panels containing data, specifically the dataTrackPanel,
     * featureTrackPanel, and headerPanel.
     */
    public void repaintDataAndHeaderPanels() {
        repaintDataAndHeaderPanels(true);
    }

    public void repaintDataPanels() {
        repaintDataAndHeaderPanels(false);
    }

    /**
     * Repaint the header and data panels.
     * <p/>
     * Note:  If running in Batch mode we force synchronous painting.  This is necessary as the
     * paint() command triggers loading of data.  If allowed to proceed asynchronously the "snapshot" batch command
     * might execute before the data from a previous command has loaded.
     *
     * @param updateCommandBar
     */
    public void repaintDataAndHeaderPanels(boolean updateCommandBar) {
        if (Globals.isBatch()) {
            Runnable r = new Runnable() {
                public void run() {
                    contentPane.revalidateDataPanels();
                    rootPane.paintImmediately(rootPane.getBounds());
                }
            };
            if (SwingUtilities.isEventDispatchThread()) {
                r.run();
            } else {
                try {
                    SwingUtilities.invokeAndWait(r);
                } catch (InterruptedException e) {
                    // Just continue
                    log.error(e);
                } catch (InvocationTargetException e) {
                    log.error(e.getMessage());
                    throw new RuntimeException(e);
                }
            }
        } else {
            contentPane.revalidateDataPanels();
            rootPane.repaint();
        }

        if (updateCommandBar) {
            contentPane.updateCurrentCoordinates();
        }
    }

    public void repaintNamePanels() {
        for (TrackPanel tp : getTrackPanels()) {
            tp.getScrollPane().getNamePanel().repaint();
        }

    }

    public void selectGenomeFromList(String genomeId) {
        contentPane.getCommandBar().selectGenome(genomeId);
    }

    public Collection<String> getSelectableGenomeIDs() {
        return contentPane.getCommandBar().getSelectableGenomeIDs();
    }


    public void doDefineGenome(ProgressMonitor monitor) {

        ProgressBar.ProgressDialog progressDialog = null;
        File archiveFile = null;

        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            GenomeBuilderDialog genomeBuilderDialog = new GenomeBuilderDialog(mainFrame, this);
            genomeBuilderDialog.setVisible(true);

            File genomeZipFile = genomeBuilderDialog.getArchiveFile();
            if (genomeBuilderDialog.isCanceled() || genomeZipFile == null) {
                return;
            }

            if (monitor != null) {
                progressDialog = ProgressBar.showProgressDialog(mainFrame, "Defining Genome...", monitor, false);
            }

            String cytobandFileName = genomeBuilderDialog.getCytobandFileName();
            String geneAnnotFileName = genomeBuilderDialog.getGeneAnnotFileName();
            String fastaFileName = genomeBuilderDialog.getFastaFileName();
            String chrAliasFile = genomeBuilderDialog.getChrAliasFileName();
            String genomeDisplayName = genomeBuilderDialog.getGenomeDisplayName();
            String genomeId = genomeBuilderDialog.getGenomeId();

            GenomeListItem genomeListItem = getGenomeManager().defineGenome(
                    genomeZipFile, cytobandFileName, geneAnnotFileName,
                    fastaFileName, chrAliasFile, genomeDisplayName,
                    genomeId, monitor);

            if (genomeListItem != null) {
                contentPane.getCommandBar().refreshGenomeListComboBox();
                contentPane.getCommandBar().selectGenome(genomeListItem.getId());
            }
            if (monitor != null) {
                monitor.fireProgressChange(100);
            }

        } catch (MaximumContigGenomeException e) {

            String genomePath = "";
            if (archiveFile != null) {
                genomePath = archiveFile.getAbsolutePath();
            }

            log.error("Failed to define genome: " + genomePath, e);

            JOptionPane.showMessageDialog(mainFrame, "Failed to define the current genome " +
                    genomePath + "\n" + e.getMessage());
        } catch (GenomeException e) {
            log.error("Failed to define genome.", e);
            MessageUtils.showMessage(e.getMessage());
        } catch (Exception e) {
            String genomePath = "";
            if (archiveFile != null) {
                genomePath = archiveFile.getAbsolutePath();
            }

            log.error("Failed to define genome: " + genomePath, e);
            MessageUtils.showMessage("Unexpected error while importing a genome: " + e.getMessage());
        } finally {
            if (progressDialog != null) {
                progressDialog.setVisible(false);
            }
            WaitCursorManager.removeWaitCursor(token);
        }
    }

    public GenomeListItem getGenomeSelectedInDropdown() {
        return contentPane.getCommandBar().getGenomeSelectedInDropdown();
    }

    /**
     * Gets the collection of genome display names currently in use.
     *
     * @return Set of display names.
     */
    public Collection<String> getGenomeDisplayNames() {
        return contentPane.getCommandBar().getGenomeDisplayNames();
    }

    void loadGenomeFromServerAction() {

        Runnable showDialog = new Runnable() {
            @Override
            public void run() {
                Collection<GenomeListItem> inputListItems = GenomeManager.getInstance().getServerGenomeArchiveList();
                if (inputListItems == null) {
                    //Not necessary to display a message, getServerGenomeArchiveList does it already
                    //IOException exc = new IOException("Unable to reach genome server");
                    //MessageUtils.showErrorMessage(exc.getMessage(), exc);
                    return;
                }

                GenomeSelectionDialog dialog = new GenomeSelectionDialog(IGV.getMainFrame(), inputListItems, ListSelectionModel.SINGLE_SELECTION);
                dialog.setVisible(true);
                List<GenomeListItem> selectedValues = dialog.getSelectedValuesList();

                if (selectedValues != null && selectedValues.size() >= 1) {

                    if (selectedValues.size() == 1 && dialog.downloadSequence()) {
                        GenomeListItem oldItem = selectedValues.get(0);
                        GenomeSelectionDialog.downloadGenome(getMainFrame(), oldItem);

                        File newLocation = new File(DirectoryManager.getGenomeCacheDirectory().getAbsolutePath(), Utilities.getFileNameFromURL(oldItem.getLocation()));

                        GenomeListItem newItem = new GenomeListItem(oldItem.getDisplayableName(), newLocation.getAbsolutePath(), oldItem.getId());
                        //Checking to see if it has a downloaded sequence might seem redundant,
                        //but if the user cancels a download we want to use the oldItem
                        if (newItem.hasDownloadedSequence()) {
                            selectedValues = Arrays.asList(newItem);
                        }
                    }

                    if (selectedValues.size() > 0) {
                        GenomeManager.getInstance().addGenomeItems(selectedValues, false);
                        getContentPane().getCommandBar().refreshGenomeListComboBox();
                        selectGenomeFromList(selectedValues.get(0).getId());
                    }
                }
            }
        };
        LongRunningTask.submit(showDialog);
    }

    /**
     * Load a .genome file directly.  This method really belongs in IGVMenuBar.
     *
     * @param monitor
     * @return
     */
    public void doLoadGenome(ProgressMonitor monitor) {

        ProgressBar.ProgressDialog progressDialog = null;
        File file = null;
        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            File importDirectory = PreferenceManager.getInstance().getLastGenomeImportDirectory();
            if (importDirectory == null) {
                PreferenceManager.getInstance().setLastGenomeImportDirectory(DirectoryManager.getUserDirectory());
            }

            // Display the dialog
            file = FileDialogUtils.chooseFile("Load Genome", importDirectory, FileDialog.LOAD);

            // If a file selection was made
            if (file != null) {
                if (monitor != null) {
                    progressDialog = ProgressBar.showProgressDialog(mainFrame, "Loading Genome...", monitor, false);
                }

                loadGenome(file.getAbsolutePath(), monitor, true);

            }
        } catch (IOException e) {
            MessageUtils.showMessage("<html>Error loading: " + file.getAbsolutePath() + "<br>" + e.getMessage());
            log.error("Error loading: " + file.getAbsolutePath(), e);
        } finally {
            WaitCursorManager.removeWaitCursor(token);
            if (monitor != null) {
                monitor.fireProgressChange(100);
            }

            if (progressDialog != null) {
                progressDialog.setVisible(false);
            }
        }

    }

    public void loadGenomeById(String genomeId) {
        if (ParsingUtils.pathExists(genomeId)) {
            try {
                IGV.getInstance().loadGenome(genomeId, null, false);
            } catch (IOException e) {
                log.error("Error loading genome file: " + genomeId, e);
            }
        } else {
            contentPane.getCommandBar().selectGenome(genomeId);
        }
    }

    public void loadGenome(String path, ProgressMonitor monitor, boolean userDefined) throws IOException {
        loadGenome(path, monitor, true, userDefined);
    }

    /**
     * @param path
     * @param monitor
     * @param addGenomeTrack Whether to display the gene track as well
     * @throws IOException
     */
    public void loadGenome(String path, ProgressMonitor monitor, boolean addGenomeTrack, boolean userDefined) throws IOException {

        File file = new File(path);
        if (file.exists()) {
            File directory = file.getAbsoluteFile().getParentFile();
            PreferenceManager.getInstance().setLastGenomeImportDirectory(directory);
        }

        resetSession(null);

        Genome genome = getGenomeManager().loadGenome(path, monitor, addGenomeTrack);
        //If genome loading cancelled
        if (genome == null) return;

        final String name = genome.getDisplayName();
        final String id = genome.getId();

        GenomeListItem genomeListItem = new GenomeListItem(name, path, id);
        getGenomeManager().addGenomeItem(genomeListItem, userDefined);

        IGVCommandBar cmdBar = contentPane.getCommandBar();
        cmdBar.refreshGenomeListComboBox();
        cmdBar.selectGenome(genomeListItem.getId());
        cmdBar.updateChromosFromGenome(genome);

        FrameManager.getDefaultFrame().setChromosomeName(genome.getHomeChromosome(), true);

        //TODO Should use EventBus/events for changing genome, clean this up a lot
        menuBar.createFileMenu();
    }


    public void enableExtrasMenu() {

        menuBar.enableExtrasMenu();
    }

    /**
     * Load a collection of tracks in a background thread.
     * <p/>
     * Note: Most of the code here is to adjust the scrollbars and split pane after loading
     *
     * @param locators
     */
    public Future loadTracks(final Collection<ResourceLocator> locators) {

        contentPane.getStatusBar().setMessage("Loading ...");

        log.debug("Run loadTracks");

        Future toRet = null;
        if (locators != null && !locators.isEmpty()) {

            // NOTE:  this work CANNOT be done on the dispatch thread, it will potentially cause deadlock if
            // dialogs are opened or other Swing tasks are done.

            NamedRunnable runnable = new NamedRunnable() {
                public void run() {

                    //Collect size statistics before loading
                    List<Map<TrackPanelScrollPane, Integer>> trackPanelAttrs = getTrackPanelAttrs();

                    loadResources(locators);

                    resetPanelHeights(trackPanelAttrs.get(0), trackPanelAttrs.get(1));

                    showLoadedTrackCount();
                }

                public String getName() {
                    return "Load Tracks";
                }
            };

            toRet = LongRunningTask.submit(runnable);
        }
        log.debug("Finish loadTracks");
        return toRet;
    }

    /**
     * Cet current track count per panel.  Needed to detect which panels
     * changed.  Also record panel sizes
     *
     * @return A 2 element list: 0th element is a map from scrollpane -> number of tracks,
     * 1st element is a map from scrollpane -> track height (in pixels)
     */
    public List<Map<TrackPanelScrollPane, Integer>> getTrackPanelAttrs() {
        Map<TrackPanelScrollPane, Integer> trackCountMap = new HashMap();
        Map<TrackPanelScrollPane, Integer> panelSizeMap = new HashMap();
        for (TrackPanel tp : getTrackPanels()) {
            TrackPanelScrollPane sp = tp.getScrollPane();
            trackCountMap.put(sp, sp.getDataPanel().getAllTracks().size());
            panelSizeMap.put(sp, sp.getDataPanel().getHeight());
        }
        return Arrays.asList(trackCountMap, panelSizeMap);
    }

    /**
     * Recalculate and set heights of track panels, based on newly loaded tracks
     *
     * @param trackCountMap scrollpane -> number of tracks
     * @param panelSizeMap  scrollpane -> height in pixels
     */
    public void resetPanelHeights(Map<TrackPanelScrollPane, Integer> trackCountMap, Map<TrackPanelScrollPane, Integer> panelSizeMap) {

        double totalHeight = 0;
        for (TrackPanel tp : getTrackPanels()) {
            TrackPanelScrollPane sp = tp.getScrollPane();
            if (trackCountMap.containsKey(sp)) {
                int prevTrackCount = trackCountMap.get(sp);
                if (prevTrackCount != sp.getDataPanel().getAllTracks().size()) {
                    int scrollPosition = panelSizeMap.get(sp);
                    if (prevTrackCount != 0 && sp.getVerticalScrollBar().isShowing()) {
                        sp.getVerticalScrollBar().setMaximum(sp.getDataPanel().getHeight());
                        sp.getVerticalScrollBar().setValue(scrollPosition);
                    }
                }
            }
            // Give a maximum "weight" of 300 pixels to each panel.  If there are no tracks, give zero
            if (sp.getTrackPanel().getTracks().size() > 0)
                totalHeight += Math.min(300, sp.getTrackPanel().getPreferredPanelHeight());
        }

        // Adjust dividers for data panel.  The data panel divider can be
        // zero if there are no data tracks loaded.
        final JideSplitPane centerSplitPane = contentPane.getMainPanel().getCenterSplitPane();
        int htotal = centerSplitPane.getHeight();
        int y = 0;
        int i = 0;
        for (Component c : centerSplitPane.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                final TrackPanel trackPanel = ((TrackPanelScrollPane) c).getTrackPanel();
                if (trackPanel.getTracks().size() > 0) {
                    int panelWeight = Math.min(300, trackPanel.getPreferredPanelHeight());
                    int dh = (int) ((panelWeight / totalHeight) * htotal);
                    y += dh;
                }
                centerSplitPane.setDividerLocation(i, y);
                i++;
            }
        }

        contentPane.getMainPanel().invalidate();
    }

    public void setGeneList(GeneList geneList) {
        setGeneList(geneList, true);
    }

    public void setGeneList(final GeneList geneList, final boolean recordHistory) {

        final CursorToken token = WaitCursorManager.showWaitCursor();

        SwingUtilities.invokeLater(new NamedRunnable() {
            public void run() {
                try {
                    if (geneList == null) {
                        session.setCurrentGeneList(null);
                    } else {
                        if (recordHistory) {
                            session.getHistory().push("List: " + geneList.getName(), 0);
                        }
                        session.setCurrentGeneList(geneList);
                    }
                    Preloader.preload();
                    resetFrames();
                } finally {
                    WaitCursorManager.removeWaitCursor(token);

                }
            }

            public String getName() {
                return "Set gene list";
            }
        });
    }

    public void setDefaultFrame(String searchString) {
        FrameManager.setToDefaultFrame(searchString);
        resetFrames();
    }

    public void resetFrames() {
        contentPane.getMainPanel().headerPanelContainer.createHeaderPanels();
        for (TrackPanel tp : getTrackPanels()) {
            tp.createDataPanels();
        }

        contentPane.getCommandBar().setGeneListMode(FrameManager.isGeneListMode());
        contentPane.getMainPanel().applicationHeaderPanel.revalidate();
        contentPane.getMainPanel().validate();
        contentPane.getMainPanel().repaint();
    }

    final public void doViewPreferences() {
        doViewPreferences(null);
    }

    /**
     * Open the user preferences dialog
     */
    final public void doViewPreferences(final String tabToSelect) {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {

                boolean originalSingleTrackValue =
                        PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY);

                PreferencesEditor dialog = new PreferencesEditor(mainFrame, true);
                if (tabToSelect != null) {
                    dialog.selectTab(tabToSelect);
                }
                dialog.setVisible(true);


                if (dialog.isCanceled()) {
                    resetStatusMessage();
                    return;

                }


                try {

                    //Should data and feature panels be combined ?
                    boolean singlePanel = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY);
                    if (originalSingleTrackValue != singlePanel) {
                        JOptionPane.showMessageDialog(mainFrame, "Panel option change will take affect after restart.");
                    }


                } finally {

                    // Update the state of the current tracks for drawing purposes
                    doRefresh();
                    resetStatusMessage();

                }


            }
        });
    }

    final public void saveStateForExit() {

        // Store recent sessions
        if (!getRecentSessionList().isEmpty()) {

            int size = getRecentSessionList().size();
            if (size > UIConstants.NUMBER_OF_RECENT_SESSIONS_TO_LIST) {
                size = UIConstants.NUMBER_OF_RECENT_SESSIONS_TO_LIST;
            }

            String recentSessions = "";
            for (int i = 0; i <
                    size; i++) {
                recentSessions += getRecentSessionList().get(i);

                if (i < (size - 1)) {
                    recentSessions += ";";
                }

            }
            PreferenceManager.getInstance().remove(PreferenceManager.RECENT_SESSION_KEY);
            PreferenceManager.getInstance().setRecentSessions(recentSessions);
        }

        // Save application location and size
        PreferenceManager.getInstance().setApplicationFrameBounds(mainFrame.getBounds());
        PreferenceManager.getInstance().put(PreferenceManager.FRAME_STATE_KEY, "" + mainFrame.getExtendedState());

    }

    final public void doShowAttributeDisplay(boolean enableAttributeView) {

        boolean oldState = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY);

        // First store the newly requested state
        PreferenceManager.getInstance().setShowAttributeView(enableAttributeView);

        //menuItem.setSelected(enableAttributeView);

        // Now, if the state has actually change we
        // need to refresh everything
        if (oldState != enableAttributeView) {
            doRefresh();
        }


    }


    final public void doRefresh() {

        contentPane.getMainPanel().revalidate();
        mainFrame.repaint();
        //getContentPane().repaint();
        contentPane.getCommandBar().updateComponentStates();
        menuBar.createFileMenu();
    }

    final public void refreshCommandBar() {
        contentPane.getCommandBar().updateCurrentCoordinates();
    }


    // TODO -- move all of this attribute stuff out of IGV,  perhaps to
    // some Attribute helper class.

    final public void doSelectDisplayableAttribute() {

        List<String> allAttributes = AttributeManager.getInstance().getAttributeNames();
        Set<String> hiddenAttributes = IGV.getInstance().getSession().getHiddenAttributes();
        final CheckListDialog dlg = new CheckListDialog(mainFrame, allAttributes, hiddenAttributes, false);
        dlg.setVisible(true);

        if (!dlg.isCanceled()) {
            // If any "default" attributes are checked turn off hide default option
            Set<String> selections = dlg.getSelections();
            for (String att : AttributeManager.defaultTrackAttributes) {
                if (selections.contains(att)) {
                    PreferenceManager.getInstance().put(PreferenceManager.SHOW_DEFAULT_TRACK_ATTRIBUTES, true);
                    break;
                }
            }
            IGV.getInstance().getSession().setHiddenAttributes(dlg.getNonSelections());
            doRefresh();
        }
    }


    final public void saveImage(Component target) {
        saveImage(target, "igv_snapshot");
    }

    final public void saveImage(Component target, String title) {
        contentPane.getStatusBar().setMessage("Creating image...");
        File defaultFile = new File(title + ".png");
        createSnapshot(target, defaultFile);
    }

    public boolean isExportingSnapshot() {
        return isExportingSnapshot;
    }

    final public void createSnapshot(final Component target, final File defaultFile) {

        File file = selectSnapshotFile(defaultFile);
        if (file == null) {
            return;
        }

        CursorToken token = null;
        try {
            token = WaitCursorManager.showWaitCursor();
            contentPane.getStatusBar().setMessage("Exporting image: " + defaultFile.getAbsolutePath());
            String msg = createSnapshotNonInteractive(target, file, false);
            if (msg != null && msg.toLowerCase().startsWith("error")) {
                MessageUtils.showMessage(msg);
            }
        } catch (IOException e) {
            log.error("Error creating exporting image ", e);
            MessageUtils.showMessage(("Error creating the image file: " + defaultFile + "<br> "
                    + e.getMessage()));
        } finally {
            if (token != null) WaitCursorManager.removeWaitCursor(token);
            resetStatusMessage();
        }

    }

    /**
     * Create a snapshot image of {@code target} and save it to {@code file}. The file type of the exported
     * snapshot will be chosen by the extension of {@code file}, which must be a supported type.
     *
     * @param target
     * @param file
     * @param paintOffscreen Whether to include offScreen data in the snapshot. Components must implement
     *                       the {@link Paintable} interface for this to work
     * @throws IOException
     * @api
     * @see SnapshotFileChooser.SnapshotFileType
     */
    public String createSnapshotNonInteractive(Component target, File file, boolean paintOffscreen) throws IOException {

        log.debug("Creating snapshot: " + file.getName());

        String extension = FileUtils.getFileExtension(file.getAbsolutePath());

        SnapshotFileChooser.SnapshotFileType type = SnapshotFileChooser.getSnapshotFileType(extension);

        String message;
        IOException exc = null;

        if (type == SnapshotFileChooser.SnapshotFileType.NULL) {
            message = "ERROR: Unknown file extension " + extension;
            log.error(message);
            return message;
        } else if (type == SnapshotFileChooser.SnapshotFileType.EPS && !SnapshotUtilities.canExportScreenshotEps()) {
            message = "ERROR: EPS output requires EPSGraphics library. See https://www.broadinstitute.org/software/igv/third_party_tools#epsgraphics";
            log.error(message);
            return message;
        }

        //boolean doubleBuffered = RepaintManager.currentManager(contentPane).isDoubleBufferingEnabled();
        try {
            setExportingSnapshot(true);
            message = SnapshotUtilities.doComponentSnapshot(target, file, type, paintOffscreen);
        } catch (IOException e) {
            exc = e;
            message = e.getMessage();
        } finally {
            setExportingSnapshot(false);
        }
        log.debug("Finished creating snapshot: " + file.getName());
        if (exc != null) throw exc;


        return message;
    }

    public File selectSnapshotFile(File defaultFile) {

        File snapshotDirectory = PreferenceManager.getInstance().getLastSnapshotDirectory();

        JFileChooser fc = new SnapshotFileChooser(snapshotDirectory, defaultFile);
        fc.showSaveDialog(mainFrame);
        File file = fc.getSelectedFile();

        // If a file selection was made
        if (file != null) {
            File directory = file.getParentFile();
            if (directory != null) {
                PreferenceManager.getInstance().setLastSnapshotDirectory(directory);
            }

        }

        return file;
    }


    private void createZoomCursors() throws HeadlessException, IndexOutOfBoundsException {
        if (zoomInCursor == null || zoomOutCursor == null) {
            final Image zoomInImage = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_IN).getImage();
            final Image zoomOutImage = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_OUT).getImage();
            final Point hotspot = new Point(10, 10);
            zoomInCursor = mainFrame.getToolkit().createCustomCursor(zoomInImage, hotspot, "Zoom in");
            zoomOutCursor = mainFrame.getToolkit().createCustomCursor(zoomOutImage, hotspot, "Zoom out");

        }

    }

    private void createHandCursor() throws HeadlessException, IndexOutOfBoundsException {
        /*if (handCursor == null) {
            BufferedImage handImage = new BufferedImage(32, 32, BufferedImage.TYPE_INT_ARGB);

            // Make backgroun transparent
            Graphics2D g = handImage.createGraphics();
            g.setComposite(AlphaComposite.getInstance(AlphaComposite.CLEAR, 0.0f));
            Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, 32, 32);
            g.fill(rect);

            // Draw hand image in middle
            g = handImage.createGraphics();
            g.drawImage(IconFactory.getInstance().getIcon(IconFactory.IconID.OPEN_HAND).getImage(), 0, 0, null);
            handCursor = getToolkit().createCustomCursor(handImage, new Point(8, 6), "Move");
        }*/

        if (fistCursor == null) {
            final BufferedImage handImage = new BufferedImage(32, 32, BufferedImage.TYPE_INT_ARGB);

            // Make backgroun transparent
            Graphics2D g = handImage.createGraphics();
            g.setComposite(AlphaComposite.getInstance(
                    AlphaComposite.CLEAR, 0.0f));
            Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, 32, 32);
            g.fill(rect);

            // Draw hand image in middle
            g = handImage.createGraphics();
            boolean ready = g.drawImage(IconFactory.getInstance().getIcon(IconFactory.IconID.FIST).getImage(), 0, 0, new ImageObserver() {
                @Override
                public boolean imageUpdate(Image img, int infoflags, int x, int y, int width, int height) {
                    if ((infoflags & ImageObserver.ALLBITS) != 0) {
                        // Image is ready
                        try {
                            fistCursor = mainFrame.getToolkit().createCustomCursor(handImage, new Point(8, 6), "Move");
                        } catch (Exception e) {
                            log.error("Could not create fistCursor", e);
                            fistCursor = Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR);
                        }
                        return false;
                    } else {
                        return true;
                    }
                }
            });
            if (ready) {
                try {
                    fistCursor = mainFrame.getToolkit().createCustomCursor(handImage, new Point(8, 6), "Move");
                } catch (Exception e) {
                    log.info("Warning: could not create fistCursor");
                    fistCursor = Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR);
                }

            }

        }

    }

    private void createDragAndDropCursor()
            throws HeadlessException, IndexOutOfBoundsException {

        if (dragNDropCursor == null) {
            ImageIcon icon = IconFactory.getInstance().getIcon(IconFactory.IconID.DRAG_AND_DROP);

            int width = icon.getIconWidth();
            int height = icon.getIconHeight();

            final BufferedImage dragNDropImage =
                    new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

            // Make background transparent
            Graphics2D g = dragNDropImage.createGraphics();
            g.setComposite(AlphaComposite.getInstance(AlphaComposite.CLEAR, 0.0f));
            Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, width, height);
            g.fill(rect);

            // Draw DND image
            g = dragNDropImage.createGraphics();
            Image image = icon.getImage();
            boolean ready = g.drawImage(image, 0, 0, new ImageObserver() {
                @Override
                public boolean imageUpdate(Image img, int infoflags, int x, int y, int width, int height) {
                    if ((infoflags & ImageObserver.ALLBITS) != 0) {
                        // Image is ready
                        try {
                            dragNDropCursor = mainFrame.getToolkit().createCustomCursor(
                                    dragNDropImage, new Point(0, 0), "Drag and Drop");
                        } catch (Exception e) {
                            log.info("Warning: could not create dragNDropCursor");
                            dragNDropCursor = Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR);
                        }
                        return false;

                    } else {
                        return true;
                    }
                }
            });
            if (ready) {
                try {
                    dragNDropCursor = mainFrame.getToolkit().createCustomCursor(
                            dragNDropImage, new Point(0, 0), "Drag and Drop");
                } catch (Exception e) {
                    log.info("Warning: could not create dragNDropCursor");
                    dragNDropCursor = Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR);
                }
            }
        }
    }

    /**
     * Set the session to the file specified by {@code sessionPath}
     * If you want to create a new session, consider {@link #newSession()}
     * as that preserves the gene track.
     *
     * @param sessionPath
     */
    public void resetSession(String sessionPath) {
        System.gc();

        AttributeManager.getInstance().clearAllAttributes();

        String tile = sessionPath == null ? UIConstants.APPLICATION_NAME : sessionPath;
        mainFrame.setTitle(tile);

        menuBar.resetSessionActions();

        AttributeManager.getInstance().clearAllAttributes();

        if (session == null) {
            session = new Session(sessionPath);
        } else {
            session.reset(sessionPath);
        }

        alignmentTrackListeners.clear();
        groupListeners.clear();

        contentPane.getMainPanel().resetPanels();

        //TODO -- this is a very blunt and dangerous way to clean up -- change to close files associated with this session
        SeekableFileStream.closeAllInstances();

        doRefresh();

        System.gc();
    }

    /**
     * Creates a new IGV session, and restores the gene track afterwards.
     * For that reason, if one wishes to keep the default gene track, this method
     * should be used, rather than resetSession
     */
    public void newSession() {
        resetSession(null);
        setGenomeTracks(GenomeManager.getInstance().getCurrentGenome().getGeneTrack());
    }

    /**
     * Set the status bar message.  If the message equals "Done." intercept
     * and reset to the default "quite" message,  currently the number of tracks
     * loaded.
     *
     * @param message
     */
    public void setStatusBarMessage(String message) {
        if (message.equals("Done.")) {
            resetStatusMessage();
        }
        contentPane.getStatusBar().setMessage(message);
    }

    /**
     * Set the status bar message.  If the message equals "Done." intercept
     * and reset to the default "quite" message,  currently the number of tracks
     * loaded.
     *
     * @param message
     */
    public void setStatusBarPosition(String message) {
        contentPane.getStatusBar().setMessage2(message);
    }

    /**
     * Resets factory settings. this is not the same as reset user defaults
     * DO NOT DELETE used when debugging
     */
    public void resetToFactorySettings() {

        try {
            PreferenceManager.getInstance().clear();
            boolean isShow = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY);
            doShowAttributeDisplay(isShow);
            Preferences prefs = Preferences.userNodeForPackage(Globals.class);
            prefs.remove(DirectoryManager.IGV_DIR_USERPREF);
            doRefresh();

        } catch (Exception e) {
            String message = "Failure while resetting preferences!";
            log.error(message, e);
            MessageUtils.showMessage(message + ": " + e.getMessage());
        }

    }

    public void setFilterMatchAll(boolean value) {
        menuBar.setFilterMatchAll(value);
    }

    public boolean isFilterMatchAll() {
        return menuBar.isFilterMatchAll();
    }

    public void setFilterShowAllTracks(boolean value) {
        menuBar.setFilterShowAllTracks(value);

    }

    public boolean isFilterShowAllTracks() {
        return menuBar.isFilterShowAllTracks();
    }

    /**
     * Add a new data panel set
     */
    public TrackPanelScrollPane addDataPanel(String name) {
        return contentPane.getMainPanel().addDataPanel(name);
    }


    /**
     * Return the panel with the given name.  This is called infrequently, and doesn't need to be fast (linear
     * search is fine).
     *
     * @param name
     * @return
     */
    public TrackPanel getTrackPanel(String name) {
        for (TrackPanel sp : getTrackPanels()) {
            if (name.equals(sp.getName())) {
                return sp;
            }
        }

        // If we get this far this is a new panel
        TrackPanelScrollPane sp = addDataPanel(name);
        return sp.getTrackPanel();
    }


    /**
     * Return an ordered list of track panels.  This method is provided primarily for storing sessions, where
     * the track panels need to be stored in order.
     */
    public List<TrackPanel> getTrackPanels() {
        return contentPane.getMainPanel().getTrackPanels();
    }


    public boolean scrollToTrack(String trackName) {
        for (TrackPanel tp : getTrackPanels()) {
            if (tp.getScrollPane().getNamePanel().scrollTo(trackName)) {
                return true;
            }
        }
        return false;
    }


    public Session getSession() {
        return session;
    }

    /**
     * Restore a session file, and optionally go to a locus.  Called upon startup and from user action.
     *
     * @param sessionFile
     * @param locus
     */
    final public void doRestoreSession(final File sessionFile, final String locus) {

        if (sessionFile.exists()) {

            doRestoreSession(sessionFile.getAbsolutePath(), locus, false);

        } else {
            String message = "Session file does not exist! : " + sessionFile.getAbsolutePath();
            log.error(message);
            MessageUtils.showMessage(message);
        }

    }

    /**
     * Load a session file, possibly asynchronously (if on the event dispatch thread).
     *
     * @param sessionPath
     * @param locus
     * @param merge
     */
    public void doRestoreSession(final String sessionPath,
                                 final String locus,
                                 final boolean merge) {

        Runnable runnable = new Runnable() {
            public void run() {
                restoreSessionSynchronous(sessionPath, locus, merge);
            }
        };
        LongRunningTask.submit(runnable);
    }

    /**
     * Load a session file in the current thread.  This should not be called from the event dispatch thread.
     *
     * @param merge
     * @param sessionPath
     * @param locus
     * @return true if successful
     */
    public boolean restoreSessionSynchronous(String sessionPath, String locus, boolean merge) {
        InputStream inputStream = null;
        try {
            if (!merge) {
                // Do this first, it closes all open SeekableFileStreams.
                resetSession(sessionPath);
            }

            setStatusBarMessage("Opening session...");
            inputStream = new BufferedInputStream(ParsingUtils.openInputStreamGZ(new ResourceLocator(sessionPath)));

            boolean isUCSC = sessionPath.endsWith(".session") || sessionPath.endsWith(".session.txt");
            boolean isIndexAware = sessionPath.endsWith(".idxsession") || sessionPath.endsWith(".idxsession.txt");
            final SessionReader sessionReader = isUCSC ?
                    new UCSCSessionReader(this) :
                    (isIndexAware ? new IndexAwareSessionReader(this) : new IGVSessionReader(this));

            sessionReader.loadSession(inputStream, session, sessionPath);

            String searchText = locus == null ? session.getLocus() : locus;

            // NOTE: Nothing to do if chr == all
            if (!FrameManager.isGeneListMode() && searchText != null &&
                    !searchText.equals(Globals.CHR_ALL) && searchText.trim().length() > 0) {
                goToLocus(searchText);
            }


            mainFrame.setTitle(UIConstants.APPLICATION_NAME + " - Session: " + sessionPath);
            System.gc();


            double[] dividerFractions = session.getDividerFractions();
            if (dividerFractions != null) {
                contentPane.getMainPanel().setDividerFractions(dividerFractions);
            }
            session.clearDividerLocations();

            //If there's a RegionNavigatorDialog, kill it.
            //this could be done through the Observer that RND uses, I suppose.  Not sure that's cleaner
            RegionNavigatorDialog.destroyInstance();

            if (!getRecentSessionList().contains(sessionPath)) {
                getRecentSessionList().addFirst(sessionPath);
            }
            doRefresh();
            return true;

        } catch (Exception e) {
            String message = "Error loading session session : <br>&nbsp;&nbsp;" + sessionPath + "<br>" +
                    e.getMessage();
            log.error(message, e);
            throw new RuntimeException(e);
        } finally {
            if (inputStream != null) {
                try {
                    inputStream.close();
                } catch (IOException iOException) {
                    log.error("Error closing session stream", iOException);
                }
                resetStatusMessage();
            }
        }
    }


    /**
     * Uses either current session.getPersistent, or preferences, depending
     * on if IGV has an instance or not. Generally intended for testing
     *
     * @param key
     * @param def
     * @return
     * @see Session#getPersistent(String, String)
     * @see PreferenceManager#getPersistent(String, String)
     */
    public static String getPersistent(String key, String def) {
        if (IGV.hasInstance()) {
            return IGV.getInstance().getSession().getPersistent(key, def);
        } else {
            return PreferenceManager.getInstance().getPersistent(key, def);
        }
    }


    /**
     * Reset the default status message, which is the number of tracks loaded.
     */
    public void resetStatusMessage() {
        contentPane.getStatusBar().setMessage("" +
                getVisibleTrackCount() + " tracks loaded");

    }

    public void rebuildGenomeDropdownList() {
        GenomeManager.getInstance().buildGenomeItemList();
        contentPane.getCommandBar().refreshGenomeListComboBox();
    }

    public void showLoadedTrackCount() {

        final int visibleTrackCount = getVisibleTrackCount();
        contentPane.getStatusBar().setMessage("" +
                visibleTrackCount + (visibleTrackCount == 1 ? " track" : " tracks"));
    }

    private void closeWindow(final ProgressBar.ProgressDialog progressDialog) {
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                progressDialog.setVisible(false);
            }
        });
    }

    /**
     * Jump to a locus synchronously. {@code locus} can be any valid search term,
     * including gene names. Genomic coordinates (e.g. "chr5:500-1000") are recommended
     * Used for port command options
     *
     * @param locus
     * @api
     */
    public void goToLocus(String locus) {
        contentPane.getCommandBar().searchByLocus(locus);
    }

    /**
     * To to multiple loci,  creating a new gene list if required.  This method is provided to support control of
     * multiple panels from a command or external program.
     *
     * @param loci
     */
    public void goToLociList(List<String> loci) {

        List<ReferenceFrame> frames = FrameManager.getFrames();
        if (frames.size() == loci.size()) {
            for (int i = 0; i < loci.size(); i++) {
                frames.get(i).jumpTo(new Locus(loci.get(i)));
            }
            repaint();
        } else {
            GeneList geneList = new GeneList("", loci, false);
            getSession().setCurrentGeneList(geneList);
            resetFrames();
        }

    }

    public void tweakPanelDivider() {
        contentPane.getMainPanel().tweakPanelDivider();
    }

    public void removeDataPanel(String name) {
        contentPane.getMainPanel().removeDataPanel(name);
    }

    public void layoutMainPanel() {
        contentPane.getMainPanel().doLayout();
    }

    public MainPanel getMainPanel() {
        return contentPane.getMainPanel();
    }

    public void setExportingSnapshot(boolean exportingSnapshot) {
        isExportingSnapshot = exportingSnapshot;
        if (isExportingSnapshot) {
            RepaintManager.currentManager(contentPane).setDoubleBufferingEnabled(false);
        } else {
            RepaintManager.currentManager(contentPane).setDoubleBufferingEnabled(true);
        }
    }

    public LinkedList<String> getRecentSessionList() {
        return recentSessionList;
    }

    public void setRecentSessionList(LinkedList<String> recentSessionList) {
        this.recentSessionList = recentSessionList;
    }

    public IGVContentPane getContentPane() {
        return contentPane;
    }

    public GenomeManager getGenomeManager() {
        return genomeManager;
    }

    JCheckBoxMenuItem showPeakMenuItem;
    PeakCommandBar peakCommandBar;

    public void addCommandBar(PeakCommandBar cb) {
        this.peakCommandBar = cb;
        contentPane.add(peakCommandBar);
        contentPane.invalidate();

        showPeakMenuItem = new JCheckBoxMenuItem("Show peaks toolbar");
        showPeakMenuItem.setSelected(true);
        showPeakMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                if (showPeakMenuItem.isSelected()) {
                    contentPane.add(peakCommandBar);
                } else {
                    contentPane.remove(peakCommandBar);
                }
            }
        });

        menuBar.getViewMenu().addSeparator();
        menuBar.getViewMenu().add(showPeakMenuItem);
    }

    public boolean isShowDetailsOnClick() {
        return contentPane != null && contentPane.getCommandBar().getDetailsBehavior() == IGVCommandBar.SHOW_DETAILS_BEHAVIOR.CLICK;
    }

    public boolean isShowDetailsOnHover() {
        return contentPane != null && contentPane.getCommandBar().getDetailsBehavior() == IGVCommandBar.SHOW_DETAILS_BEHAVIOR.HOVER;
    }

    public void openStatusWindow() {
        if (statusWindow == null) {
            statusWindow = new StatusWindow();
        }
        statusWindow.setVisible(true);
    }

    public void setStatusWindowText(String text) {
        if (statusWindow != null && statusWindow.isVisible()) {
            statusWindow.updateText(text);
        }
    }


    /**
     * Load resources into IGV. Tracks are added to the appropriate panel
     *
     * @param locators
     */
    public void loadResources(Collection<ResourceLocator> locators) {

        //Set<TrackPanel> changedPanels = new HashSet();

        log.info("Loading " + locators.size() + " resources.");
        final MessageCollection messages = new MessageCollection();


        // Load files concurrently -- TODO, put a limit on # of threads?
        List<Thread> threads = new ArrayList<Thread>(locators.size());

        for (final ResourceLocator locator : locators) {

            // If its a local file, check explicitly for existence (rather than rely on exception)
            if (locator.isLocal()) {
                File trackSetFile = new File(locator.getPath());
                if (!trackSetFile.exists()) {
                    messages.append("File not found: " + locator.getPath() + "\n");
                    continue;
                }
            }

            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        List<Track> tracks = load(locator);
                        log.debug(tracks.size() + " new tracks loaded");
                        addTracks(tracks, locator);
                    } catch (Exception e) {
                        log.error("Error loading track", e);
                        messages.append("Error loading " + locator + ": " + e.getMessage());
                    }
                }
            };

            //Thread thread = new Thread(runnable);
            //thread.start();
            //threads.add(thread);
            runnable.run();
        }

        // Wait for all threads to complete
        for (Thread t : threads) {
            try {
                t.join();
            } catch (InterruptedException ignore) {
                log.error(ignore.getMessage(), ignore);
                messages.append("Thread interrupted: " + ignore.getMessage());
            }
        }

        resetGroups();
        resetOverlayTracks();

        if (!messages.isEmpty()) {
            for (String message : messages.getMessages()) {
                MessageUtils.showMessage(message);
            }
        }
    }

    /**
     * Add tracks to the specified panel
     *
     * @param tracks
     * @param panelName
     * @api
     */
    public void addTracks(List<Track> tracks, PanelName panelName) {
        TrackPanel panel = getTrackPanel(panelName.getName());
        panel.addTracks(tracks);
        doRefresh();
    }

    /**
     * Add the specified tracks to the appropriate panel. Panel
     * is chosen based on characteristics of the {@code locator}.
     *
     * @param tracks
     * @param locator
     */
    void addTracks(List<Track> tracks, ResourceLocator locator) {
        if (tracks.size() > 0) {
            String path = locator.getPath();

            // Get an appropriate panel.  If its a VCF file create a new panel if the number of genotypes
            // is greater than 10
            TrackPanel panel = getPanelFor(locator);
            if (path.endsWith(".vcf") || path.endsWith(".vcf.gz") ||
                    path.endsWith(".vcf4") || path.endsWith(".vcf4.gz")) {
                Track t = tracks.get(0);
                if (t instanceof VariantTrack && ((VariantTrack) t).getAllSamples().size() > 10) {
                    String newPanelName = "Panel" + System.currentTimeMillis();
                    panel = addDataPanel(newPanelName).getTrackPanel();
                }
            }
            panel.addTracks(tracks);
        }
    }


    /**
     * Load a resource and return the tracks.
     * Does not automatically add anything
     *
     * @param locator
     * @return A list of loaded tracks
     */
    public List<Track> load(ResourceLocator locator) throws DataLoadException {

        TrackLoader loader = new TrackLoader();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        List<Track> newTracks = loader.load(locator, genome);
        if (newTracks.size() > 0) {
            for (Track track : newTracks) {
                String fn = locator.getPath();
                int lastSlashIdx = fn.lastIndexOf("/");
                if (lastSlashIdx < 0) {
                    lastSlashIdx = fn.lastIndexOf("\\");
                }
                if (lastSlashIdx > 0) {
                    fn = fn.substring(lastSlashIdx + 1);
                }
                track.setAttributeValue(Globals.TRACK_NAME_ATTRIBUTE, track.getName());
                track.setAttributeValue(Globals.TRACK_DATA_FILE_ATTRIBUTE, fn);
                track.setAttributeValue(Globals.TRACK_DATA_TYPE_ATTRIBUTE, track.getTrackType().toString());

                // If the track listens for group events add it to the listener list
                if (track instanceof TrackGroupEventListener) {
                    addGroupEventListener((TrackGroupEventListener) track);
                }
            }
        }

        return newTracks;
    }


    /**
     * Load the data file into the specified panel.   Triggered via drag and drop.
     */
    public void load(ResourceLocator locator, TrackPanel panel) throws DataLoadException {
        // If this is a session  TODO -- need better "is a session?" test
        if (locator.getPath().endsWith(".xml") || locator.getPath().endsWith(("session"))) {
            boolean merge = false;  // TODO -- ask user?
            this.doRestoreSession(locator.getPath(), null, merge);
        }

        // Not a session, load into target panel
        List<Track> tracks = load(locator);
        panel.addTracks(tracks);
        doRefresh();
    }

    /**
     * Return a DataPanel appropriate for the resource type
     *
     * @param locator
     * @return
     */
    public TrackPanel getPanelFor(ResourceLocator locator) {
        String path = locator.getPath().toLowerCase();
        if ("alist".equals(locator.getType())) {
            return getVcfBamPanel();
        } else if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY)) {
            return getTrackPanel(DATA_PANEL_NAME);
        } else if (path.endsWith(".sam") || path.endsWith(".bam") ||
                path.endsWith(".sam.list") || path.endsWith(".bam.list") ||
                path.endsWith(".aligned") || "ga4gh".equals(locator.getType())) {

            String newPanelName = "Panel" + System.currentTimeMillis();
            return addDataPanel(newPanelName).getTrackPanel();
        } else {
            return getDefaultPanel(locator);
        }
    }

    public Set<TrackType> getLoadedTypes() {
        Set<TrackType> types = new HashSet();
        for (Track t : getAllTracks()) {
            TrackType type = t.getTrackType();
            if (t != null) {
                types.add(type);
            }
        }
        return types;
    }


    /**
     * Experimental method to support VCF -> BAM coupling
     *
     * @return
     */
    public TrackPanel getVcfBamPanel() {
        String panelName = "VCF_BAM";
        TrackPanel panel = getTrackPanel(panelName);
        if (panel != null) {
            return panel;
        } else {
            return addDataPanel(panelName).getTrackPanel();
        }
    }


    private TrackPanel getDefaultPanel(ResourceLocator locator) {

        if (locator.getType() != null && locator.getType().equalsIgnoreCase("das")) {
            return getTrackPanel(FEATURE_PANEL_NAME);
        }

        String filename = locator.getPath().toLowerCase();

        if (filename.endsWith(".txt") || filename.endsWith(".tab") || filename.endsWith(
                ".xls") || filename.endsWith(".gz")) {
            filename = filename.substring(0, filename.lastIndexOf("."));
        }


        if (filename.contains("refflat") || filename.contains("ucscgene") ||
                filename.contains("genepred") || filename.contains("ensgene") ||
                filename.contains("refgene") ||
                filename.endsWith("gff") || filename.endsWith("gtf") ||
                filename.endsWith("gff3") || filename.endsWith("embl") ||
                filename.endsWith("bed") || filename.endsWith("gistic") ||
                filename.endsWith("bedz") || filename.endsWith("repmask") ||
                filename.contains("dranger")) {
            return getTrackPanel(FEATURE_PANEL_NAME);
        } else {
            return getTrackPanel(DATA_PANEL_NAME);
        }
    }

    public void reset() {
        groupByAttribute = null;
        for (TrackPanel sp : getTrackPanels()) {
            if (DATA_PANEL_NAME.equals(sp.getName())) {
                sp.reset();
                break;
            }
        }
        groupListeners.clear();
    }


    public boolean sortAlignmentTracks(AlignmentTrack.SortOption option, String tag) {
        return sortAlignmentTracks(option, null, tag);
    }

    public boolean sortAlignmentTracks(AlignmentTrack.SortOption option, Double location, String tag) {
        double actloc;
        boolean toRet = true;
        for (Track t : getAllTracks()) {
            if (t instanceof AlignmentTrack) {
                for (ReferenceFrame frame : FrameManager.getFrames()) {
                    actloc = location != null ? location : frame.getCenter();
                    toRet &= ((AlignmentTrack) t).sortRows(option, frame, actloc, tag);
                }
            }
        }
        return toRet;
    }

    /**
     * Group all alignment tracks by the specified option.
     *
     * @param option
     * @api
     */
    public void groupAlignmentTracks(AlignmentTrack.GroupOption option, String tag) {
        for (Track t : getAllTracks()) {
            if (t instanceof AlignmentTrack) {
                ((AlignmentTrack) t).groupAlignments(option, tag);
            }
        }
    }

    public void packAlignmentTracks() {
        for (Track t : getAllTracks()) {
            if (t instanceof AlignmentTrack) {
                ((AlignmentTrack) t).packAlignments();
            }
        }
    }


    public void setTrackDisplayMode(Track.DisplayMode mode, String trackName) {
        for (Track t : getAllTracks()) {
            if (trackName == null || t.getName().equals(trackName)) {
                t.setDisplayMode(mode);
            }
        }

    }


    /**
     * Reset the overlay tracks collection.  Currently the only overlayable track
     * type is Mutation.  This method finds all mutation tracks and builds a map
     * of key -> mutation track,  where the key is the specified attribute value
     * for linking tracks for overlay.
     */
    public void resetOverlayTracks() {
        log.debug("Resetting Overlay Tracks");
        overlayTracksMap.clear();
        overlaidTracks.clear();


        // Old option to allow overlaying based on an arbitrary attribute.
        // String overlayAttribute = igv.getSession().getOverlayAttribute();

        for (Track track : getAllTracks()) {
            if (track != null && track.getTrackType() == TrackType.MUTATION) {

                String sample = track.getSample();

                if (sample != null) {
                    List<Track> trackList = overlayTracksMap.get(sample);

                    if (trackList == null) {
                        trackList = new ArrayList();
                        overlayTracksMap.put(sample, trackList);
                    }

                    trackList.add(track);
                }
            }

        }

        for (Track track : getAllTracks()) {
            if (track != null) {  // <= this should not be neccessary
                if (track.getTrackType() != TrackType.MUTATION) {
                    String sample = track.getSample();
                    if (sample != null) {
                        List<Track> trackList = overlayTracksMap.get(sample);
                        if (trackList != null) overlaidTracks.addAll(trackList);
                    }
                }
            }
        }

        boolean displayOverlays = getSession().getOverlayMutationTracks();
        for (Track track : getAllTracks()) {
            if (track != null) {
                if (track.getTrackType() == TrackType.MUTATION) {
                    track.setOverlayed(displayOverlays && overlaidTracks.contains(track));
                }
            }
        }
    }


    /**
     * Return tracks overlaid on "track"
     * // TODO -- why aren't overlaid tracks stored in a track member?  This seems unnecessarily complex
     *
     * @param track
     * @return
     */
    public List<Track> getOverlayTracks(Track track) {
        String sample = track.getSample();
        if (sample != null) {
            return overlayTracksMap.get(sample);
        }
        return null;
    }

    public int getVisibleTrackCount() {
        int count = 0;
        for (TrackPanel tsv : getTrackPanels()) {
            count += tsv.getVisibleTrackCount();

        }
        return count;
    }

    /**
     * Return the list of all tracks in the order they appear on the screen
     *
     * @return
     */
    public List<Track> getAllTracks() {
        List<Track> allTracks = new ArrayList<Track>();
        for (TrackPanel tp : getTrackPanels()) {
            allTracks.addAll(tp.getTracks());
        }
        return allTracks;
    }

    public List<FeatureTrack> getFeatureTracks() {
        Iterable<FeatureTrack> featureTracksIter = Iterables.filter(getAllTracks(), FeatureTrack.class);
        List<FeatureTrack> featureTracks = Lists.newArrayList(featureTracksIter);
        return featureTracks;
    }

    public List<DataTrack> getDataTracks() {
        Iterable<DataTrack> dataTracksIter = Iterables.filter(getAllTracks(), DataTrack.class);
        List<DataTrack> dataTracks = Lists.newArrayList(dataTracksIter);
        return dataTracks;
    }

    public void clearSelections() {
        for (Track t : getAllTracks()) {
            if (t != null)
                t.setSelected(false);

        }

    }

    public void setTrackSelections(Iterable<Track> selectedTracks) {
        for (Track t : selectedTracks) {
            t.setSelected(true);
        }
    }

    public void shiftSelectTracks(Track track) {
        List<Track> allTracks = getAllTracks();
        int clickedTrackIndex = allTracks.indexOf(track);
        // Find another track that is already selected.  The semantics of this
        // are not well defined, so any track will do
        int otherIndex = clickedTrackIndex;
        for (int i = 0; i < allTracks.size(); i++) {
            if (allTracks.get(i).isSelected() && i != clickedTrackIndex) {
                otherIndex = i;
                break;
            }
        }

        int left = Math.min(otherIndex, clickedTrackIndex);
        int right = Math.max(otherIndex, clickedTrackIndex);
        for (int i = left; i <= right; i++) {
            allTracks.get(i).setSelected(true);
        }
    }

    public void toggleTrackSelections(Iterable<Track> selectedTracks) {
        for (Track t : selectedTracks) {
            t.setSelected(!t.isSelected());
        }
    }

    public List<Track> getSelectedTracks() {
        ArrayList<Track> selectedTracks = new ArrayList();
        for (Track t : getAllTracks()) {
            if (t != null && t.isSelected()) {
                selectedTracks.add(t);
            }
        }
        return selectedTracks;

    }

    /**
     * Return the complete set of unique DataResourceLocators currently loaded
     *
     * @return
     */
    public Set<ResourceLocator> getDataResourceLocators() {
        HashSet<ResourceLocator> locators = new HashSet();

        for (Track track : getAllTracks()) {
            Collection<ResourceLocator> tlocators = track.getResourceLocators();

            if (tlocators != null) {
                locators.addAll(tlocators);
            }
        }
        locators.remove(null);
        return locators;

    }


    public void setAllTrackHeights(int newHeight) {
        for (Track track : getAllTracks()) {
            track.setHeight(newHeight, true);
        }

    }


    public void removeTracks(Collection<? extends Track> tracksToRemove) {

        // Make copy of list as we will be modifying the original in the loop
        List<TrackPanel> panels = getTrackPanels();
        for (TrackPanel trackPanel : panels) {
            trackPanel.removeTracks(tracksToRemove);

            if (!trackPanel.hasTracks()) {
                removeDataPanel(trackPanel.getName());
            }
        }

        for (Track t : tracksToRemove) {
            if (t instanceof TrackGroupEventListener) {
                removeGroupEventListener((TrackGroupEventListener) t);
            }
            if (t instanceof AlignmentTrackEventListener) {
                removeAlignmentTrackEvent((AlignmentTrackEventListener) t);
            }
        }

        for (Track t : tracksToRemove) {
            t.dispose();
        }
    }

    /**
     * Add gene and sequence tracks.  This is called upon switching genomes.
     *
     * @param newGeneTrack
     * @param
     */
    public void setGenomeTracks(FeatureTrack newGeneTrack) {

        TrackPanel panel = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY) ?
                getTrackPanel(DATA_PANEL_NAME) : getTrackPanel(FEATURE_PANEL_NAME);
        SequenceTrack newSeqTrack = new SequenceTrack("Reference sequence");

        panel.addTrack(newSeqTrack);
        if (newGeneTrack != null) {
            panel.addTrack(newGeneTrack);
        }

    }

    public boolean hasGeneTrack() {
        FeatureTrack geneTrack = GenomeManager.getInstance().getCurrentGenome().getGeneTrack();
        if (geneTrack == null) return false;
        for (Track t : getFeatureTracks()) {
            if (geneTrack == t) return true;
        }
        return false;
    }

    public boolean hasSequenceTrack() {
        return getSequenceTrack() != null;
    }

    /**
     * @return First SequenceTrack found, or null if none
     */
    public SequenceTrack getSequenceTrack() {
        for (Track t : getAllTracks()) {
            if (t instanceof SequenceTrack) return (SequenceTrack) t;
        }
        return null;
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    // Sorting


    /**
     * Sort all groups (data and feature) by attribute value(s).  Tracks are
     * sorted within groups.
     *
     * @param attributeNames
     * @param ascending
     */
    public void sortAllTracksByAttributes(final String attributeNames[], final boolean[] ascending) {
        assert attributeNames.length == ascending.length;

        for (TrackPanel trackPanel : getTrackPanels()) {
            trackPanel.sortTracksByAttributes(attributeNames, ascending);
        }
    }


    /**
     * Sort all groups (data and feature) by a computed score over a region.  The
     * sort is done twice (1) groups are sorted with the featureGroup, and (2) the
     * groups themselves are sorted.
     *
     * @param region
     * @param type
     * @param frame
     */
    public void sortByRegionScore(RegionOfInterest region,
                                  final RegionScoreType type,
                                  final ReferenceFrame frame) {

        final RegionOfInterest r = region == null ? new RegionOfInterest(frame.getChrName(), (int) frame.getOrigin(),
                (int) frame.getEnd() + 1, frame.getName()) : region;

        // Create a rank order of samples.  This is done globally so sorting is consistent across groups and panels.
        final List<String> sortedSamples = sortSamplesByRegionScore(r, type, frame);

        for (TrackPanel trackPanel : getTrackPanels()) {
            trackPanel.sortByRegionsScore(r, type, frame, sortedSamples);
        }
        repaintDataPanels();
    }


    /**
     * Sort a collection of tracks by a score over a region.
     *
     * @param region
     * @param type
     * @param frame
     */
    private List<String> sortSamplesByRegionScore(final RegionOfInterest region,
                                                  final RegionScoreType type,
                                                  final ReferenceFrame frame) {

        // Get the sortable tracks for this score (data) type
        final List<Track> allTracks = getAllTracks();
        final List<Track> tracksWithScore = new ArrayList(allTracks.size());
        for (Track t : allTracks) {
            if (t.isRegionScoreType(type)) {
                tracksWithScore.add(t);
            }
        }

        // Sort the "sortable" tracks
        sortByRegionScore(tracksWithScore, region, type, frame);

        // Now get sample order from sorted tracks, use to sort (tracks which do not implement the selected "sort by" score)
        List<String> sortedSamples = new ArrayList(tracksWithScore.size());
        for (Track t : tracksWithScore) {
            String att = t.getSample(); //t.getAttributeValue(linkingAtt);
            if (att != null) {
                sortedSamples.add(att);
            }

        }

        return sortedSamples;
    }

    static void sortByRegionScore(List<Track> tracks,
                                  final RegionOfInterest region,
                                  final RegionScoreType type,
                                  ReferenceFrame frame) {
        if ((tracks != null) && (region != null) && !tracks.isEmpty()) {
            final String frameName = frame != null ? frame.getName() : null;
            int tmpzoom = frame != null ? frame.getZoom() : 0;
            final int zoom = Math.max(0, tmpzoom);
            final String chr = region.getChr();
            final int start = region.getStart();
            final int end = region.getEnd();

            Comparator<Track> c = new Comparator<Track>() {

                public int compare(Track t1, Track t2) {
                    try {
                        if (t1 == null && t2 == null) return 0;
                        if (t1 == null) return 1;
                        if (t2 == null) return -1;


                        float s1 = t1.getRegionScore(chr, start, end, zoom, type, frameName);
                        float s2 = t2.getRegionScore(chr, start, end, zoom, type, frameName);

                        return Float.compare(s2, s1);


                    } catch (Exception e) {
                        log.error("Error sorting tracks. Sort might not be accurate.", e);
                        return 0;
                    }

                }
            };
            Collections.sort(tracks, c);

        }
    }


    ////////////////////////////////////////////////////////////////////////////////////////
    // Groups

    public String getGroupByAttribute() {
        return groupByAttribute;
    }


    public void setGroupByAttribute(String attributeName) {
        groupByAttribute = attributeName;
        resetGroups();
        // Some tracks need to respond to changes in grouping, fire notification event
        notifyGroupEvent();
    }


    private void resetGroups() {
        log.debug("Resetting Groups");
        for (TrackPanel trackPanel : getTrackPanels()) {
            trackPanel.groupTracksByAttribute(groupByAttribute);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Events

    public synchronized void addGroupEventListener(TrackGroupEventListener l) {
        groupListeners.add(new SoftReference<TrackGroupEventListener>(l));
    }

    public synchronized void removeGroupEventListener(TrackGroupEventListener l) {

        for (Iterator<SoftReference<TrackGroupEventListener>> it = groupListeners.iterator(); it.hasNext(); ) {
            TrackGroupEventListener listener = it.next().get();
            if (listener != null && listener == l) {
                it.remove();
                break;
            }
        }
    }

    public void notifyGroupEvent() {
        TrackGroupEvent e = new TrackGroupEvent(this);
        for (SoftReference<TrackGroupEventListener> ref : groupListeners) {
            TrackGroupEventListener l = ref.get();
            l.onTrackGroupEvent(e);
        }
    }

    public synchronized void addAlignmentTrackEventListener(AlignmentTrackEventListener l) {
        alignmentTrackListeners.add(new SoftReference<AlignmentTrackEventListener>(l));
    }

    public synchronized void removeAlignmentTrackEvent(AlignmentTrackEventListener l) {
        for (Iterator<SoftReference<AlignmentTrackEventListener>> it = alignmentTrackListeners.iterator(); it.hasNext(); ) {
            AlignmentTrackEventListener listener = it.next().get();
            if (listener != null && listener == l) {
                it.remove();
                break;
            }
        }
    }

    public void notifyAlignmentTrackEvent(Object source, AlignmentTrackEvent.Type type, boolean value) {
        AlignmentTrackEvent e = new AlignmentTrackEvent(source, type, value);
        for (SoftReference<AlignmentTrackEventListener> ref : alignmentTrackListeners) {
            AlignmentTrackEventListener l = ref.get();
            l.onAlignmentTrackEvent(e);
        }
    }


    public void notifyAlignmentTrackEvent(Object source, AlignmentTrackEvent.Type type) {
        AlignmentTrackEvent e = new AlignmentTrackEvent(source, type);
        for (SoftReference<AlignmentTrackEventListener> ref : alignmentTrackListeners) {
            AlignmentTrackEventListener l = ref.get();
            l.onAlignmentTrackEvent(e);
        }
    }


    //////////////////////////////////////////////////////////////////////////////////////////
    // Startup


    public Future startUp(Main.IGVArgs igvArgs) {

        if (log.isDebugEnabled()) {
            log.debug("startUp");
        }

        return LongRunningTask.submit(new StartupRunnable(igvArgs));
    }

    public void setRulerEnabled(boolean rulerEnabled) {
        this.rulerEnabled = rulerEnabled;
    }

    public boolean isRulerEnabled() {
        return rulerEnabled;
    }

    /**
     * Swing worker class to startup IGV
     */
    public class StartupRunnable implements Runnable {
        Main.IGVArgs igvArgs;

        StartupRunnable(Main.IGVArgs args) {
            this.igvArgs = args;

        }

        @Override
        public void run() {

            final boolean runningBatch = igvArgs.getBatchFile() != null;
            BatchRunner.setIsBatchMode(runningBatch);

            final ProgressMonitor monitor = new ProgressMonitor();
            final ProgressBar.ProgressDialog progressDialog = ProgressBar.showProgressDialog(mainFrame, "Initializing...", monitor, false);
            progressDialog.getProgressBar().setIndeterminate(true);
            monitor.fireProgressChange(20);

            mainFrame.setIconImage(getIconImage());
            if (Globals.IS_MAC) {
                setAppleDockIcon();
            }

            final PreferenceManager preferenceManager = PreferenceManager.getInstance();

            try {
                contentPane.getCommandBar().initializeGenomeList(monitor);
            } catch (FileNotFoundException ex) {
                JOptionPane.showMessageDialog(mainFrame, "Error initializing genome list: " + ex.getMessage());
                log.error("Error initializing genome list: ", ex);
            } catch (NoRouteToHostException ex) {
                JOptionPane.showMessageDialog(mainFrame, "Network error initializing genome list: " + ex.getMessage());
                log.error("Network error initializing genome list: ", ex);
            } finally {
                monitor.fireProgressChange(50);
                closeWindow(progressDialog);
            }

            if (igvArgs.getGenomeId() != null) {
                IGV.getInstance().loadGenomeById(igvArgs.getGenomeId());
            } else if (igvArgs.getSessionFile() == null) {
                String genomeId = preferenceManager.getDefaultGenome();
                contentPane.getCommandBar().selectGenome(genomeId);
            }

            //Load plugins
            //Do this before loading files, hooks might need to be inserted
            initIGVPlugins();

            //If there is an argument assume it is a session file or url
            if (igvArgs.getSessionFile() != null || igvArgs.getDataFileString() != null) {

                if (log.isDebugEnabled()) {
                    log.debug("Loading session data");
                }

                final IndefiniteProgressMonitor indefMonitor = new IndefiniteProgressMonitor();
                final ProgressBar.ProgressDialog progressDialog2 = ProgressBar.showProgressDialog(mainFrame, "Loading session data", indefMonitor, false);
                indefMonitor.start();


                if (log.isDebugEnabled()) {
                    log.debug("Calling restore session");
                }


                if (igvArgs.getSessionFile() != null) {
                    boolean success = false;
                    if (HttpUtils.isRemoteURL(igvArgs.getSessionFile())) {
                        boolean merge = false;
                        success = restoreSessionSynchronous(igvArgs.getSessionFile(), igvArgs.getLocusString(), merge);
                    } else {
                        File sf = new File(igvArgs.getSessionFile());
                        if (sf.exists()) {
                            success = restoreSessionSynchronous(sf.getAbsolutePath(), igvArgs.getLocusString(), false);
                        }
                    }
                    if (!success) {
                        String genomeId = preferenceManager.getDefaultGenome();
                        contentPane.getCommandBar().selectGenome(genomeId);

                    }
                } else if (igvArgs.getDataFileString() != null) {
                    // Not an xml file, assume its a list of data files
                    String[] dataFiles = igvArgs.getDataFileString().split(",");
                    String[] names = null;
                    if (igvArgs.getName() != null) {
                        names = igvArgs.getName().split(",");
                    }
                    String[] indexFiles = null;
                    if (igvArgs.getIndexFile() != null) {
                        indexFiles = igvArgs.getIndexFile().split(",");
                    }
                    String[] coverageFiles = null;
                    if (igvArgs.getCoverageFile() != null) {
                        coverageFiles = igvArgs.getCoverageFile().split(",");
                    }


                    List<ResourceLocator> locators = new ArrayList();


                    for (int i = 0; i < dataFiles.length; i++) {

                        String p = dataFiles[i].trim();

                        // Decode local file paths
                        if (HttpUtils.isURL(p) && !FileUtils.isRemote(p)) {
                            p = StringUtils.decodeURL(p);
                        }

                        ResourceLocator rl = new ResourceLocator(p);

                        if (names != null && i < names.length) {
                            String name = names[i];
                            rl.setName(name);
                        }

                        //Set index file, iff one was passed
                        if (indexFiles != null && i < indexFiles.length) {
                            String idxP = indexFiles[i];
                            if (HttpUtils.isURL(idxP) && !FileUtils.isRemote(idxP)) {
                                idxP = StringUtils.decodeURL(idxP);
                            }
                            if (idxP.length() > 0) {
                                rl.setIndexPath(idxP);
                            }
                        }

                        //Set coverage file, iff one was passed
                        if (coverageFiles != null && i < coverageFiles.length) {
                            String covP = coverageFiles[i];
                            if (HttpUtils.isURL(covP) && !FileUtils.isRemote(covP)) {
                                covP = StringUtils.decodeURL(covP);
                            }
                            if (covP.length() > 0) {
                                rl.setCoverage(covP);
                            }
                        }

                        locators.add(rl);
                    }
                    loadTracks(locators);
                }


                indefMonitor.stop();
                closeWindow(progressDialog2);

            }

            session.recordHistory();

            // Start up a port listener.  Port # can be overriden with "-p" command line switch
            boolean portEnabled = preferenceManager.getAsBoolean(PreferenceManager.PORT_ENABLED);
            String portString = igvArgs.getPort();
            if (portEnabled || portString != null) {
                // Command listener thread
                int port = preferenceManager.getAsInt(PreferenceManager.PORT_NUMBER);
                if (portString != null) {
                    port = Integer.parseInt(portString);
                }
                CommandListener.start(port);
            }

            UIUtilities.invokeOnEventThread(new Runnable() {
                public void run() {
                    mainFrame.setVisible(true);
                    if (igvArgs.getLocusString() != null) {
                        goToLocus(igvArgs.getLocusString());
                    }
                    if (runningBatch) {
                        LongRunningTask.submit(new BatchRunner(igvArgs.getBatchFile()));
                    } else {
//                        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_LOS)) {
//                            (new LOSDialog(mainFrame)).setVisible(true);
//                        }
                    }

                }
            });

            synchronized (IGV.getInstance()) {
                IGV.getInstance().notifyAll();
            }
        }

        private void setAppleDockIcon() {
            try {
                Image image = getIconImage();
                OSXAdapter.setDockIconImage(image);
            } catch (Exception e) {
                log.error("Error setting apple dock icon", e);
            }
        }

        private Image getIconImage() {
            String path = "resources/IGV_64.png";
            URL url = IGV.class.getResource(path);
            Image image = new ImageIcon(url).getImage();
            return image;
        }


        private void initIGVPlugins() {
            List<String> pluginClassNames = new ArrayList<String>(2);
            InputStream is = IGV.class.getResourceAsStream("resources/builtin_plugin_list.txt");
            if (is != null) {
                BufferedReader br = new BufferedReader(new InputStreamReader(is));
                String line = null;
                try {
                    while ((line = br.readLine()) != null) {
                        if (line.startsWith("##")) continue;
                        pluginClassNames.add(line);
                    }
                } catch (IOException e) {
                    log.error("Error reading builtin plugin list", e);
                }
            }
            pluginClassNames.addAll(Arrays.asList(PreferenceManager.getInstance().getIGVPluginList()));
            for (String classname : pluginClassNames) {
                if (classname.startsWith("#")) continue;
                try {
                    Class clazz = Class.forName(classname);
                    IGVPlugin plugin = (IGVPlugin) clazz.newInstance();
                    plugin.init();
                } catch (Exception e) {
                    log.error("Error loading " + classname, e);
                }
            }

        }

    }

    public static void copySequenceToClipboard(Genome genome, String chr, int start, int end) {
        try {
            IGV.getMainFrame().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
            byte[] seqBytes = genome.getSequence(chr, start, end);

            if (seqBytes == null) {
                MessageUtils.showMessage("Sequence not available. Try enabling http byte-range requests");
            } else {
                String sequence = new String(seqBytes);
                //TODO This will complement sequence if sequence track is flipped
                //Might be un-intuitive to user if they do it from region dialog
//                SequenceTrack sequenceTrack = IGV.getInstance().getSequenceTrack();
//                if(sequenceTrack != null && sequenceTrack.getStrand() == Strand.NEGATIVE){
//                    sequence = AminoAcidManager.getNucleotideComplement(sequence);
//                }
                StringUtils.copyTextToClipboard(sequence);
            }

        } finally {
            IGV.getMainFrame().setCursor(Cursor.getDefaultCursor());
        }
    }


    /**
     * Wrapper for igv.wait(timeout)
     *
     * @param timeout
     * @return True if method completed before interruption (not necessarily before timeout), otherwise false
     */
    public boolean waitForNotify(long timeout) {
        boolean completed = false;
        synchronized (this) {
            while (!completed) {
                try {
                    this.wait(timeout);
                    completed = true;
                } catch (InterruptedException e) {

                }
                break;
            }
        }
        return completed;
    }


}
