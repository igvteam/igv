/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2018 Broad Institute
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

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.jidesoft.swing.JideSplitPane;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.annotations.ForTesting;
import org.broad.igv.batch.BatchRunner;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.event.*;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.*;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.*;
import org.broad.igv.ga4gh.OAuthUtils;
import org.broad.igv.lists.GeneList;
import org.broad.igv.peaks.PeakCommandBar;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferenceEditorFX;
import org.broad.igv.prefs.PreferencesEditor;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.InsertionSelectionEvent;
import org.broad.igv.session.*;
import org.broad.igv.track.*;
import org.broad.igv.ui.WaitCursorManager.CursorToken;
import org.broad.igv.ui.dnd.GhostGlassPane;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.*;
import org.broad.igv.util.*;
import org.broad.igv.variant.VariantTrack;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.List;
import java.util.concurrent.Future;
import java.util.prefs.Preferences;

import static org.broad.igv.prefs.Constants.*;

/**
 * Represents an IGV instance, consisting of a main window and associated model.
 *
 * @author jrobinso
 */
public class IGV implements IGVEventObserver {

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

    private List<JComponent> otherToolMenus = new ArrayList<>();

    // Vertical line that follows the mouse
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

        final IGVPreferences preferences = PreferencesManager.getPreferences();

        genomeManager = GenomeManager.getInstance();
        mainFrame = frame;
        mainFrame.addWindowListener(new WindowAdapter() {


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
       // consumeEvents(glassPane);

        dNdGlassPane = new GhostGlassPane();

        mainFrame.pack();

        //Certain components MUST be visible, so we set minimum size
        //{@link MainPanel#addDataPanel}
        mainFrame.setMinimumSize(new Dimension(300, 300));

        // Set the application's previous location and size
        Dimension screenBounds = Toolkit.getDefaultToolkit().getScreenSize();
        Rectangle applicationBounds = preferences.getApplicationFrameBounds();

        if (applicationBounds == null || applicationBounds.getMaxX() > screenBounds.getWidth() ||
                applicationBounds.getMaxY() > screenBounds.getHeight() ||
                applicationBounds.width == 0 || applicationBounds.height == 0) {
            int width = Math.min(1150, (int) screenBounds.getWidth());
            int height = Math.min(800, (int) screenBounds.getHeight());
            applicationBounds = new Rectangle(0, 0, width, height);
        }
        mainFrame.setBounds(applicationBounds);

        subscribeToEvents();
    }

    private void consumeEvents(Component glassPane) {

        glassPane.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {

                Point glassPanePoint = e.getPoint();
                Container container = IGV.this.contentPane;
                Point containerPoint = SwingUtilities.convertPoint(glassPane,
                        glassPanePoint, container);

                Component component = SwingUtilities.getDeepestComponentAt(
                        container, containerPoint.x, containerPoint.y);

                if (component == IGV.this.contentPane.getStatusBar().stopButton) {
                    IGVEventBus.getInstance().post(new StopEvent());
                }
                e.consume();

            }

            @Override
            public void mousePressed(MouseEvent e) {
                e.consume();

            }
        });
        glassPane.setFocusable(true);
        glassPane.addKeyListener(new KeyListener() {
            @Override
            public void keyTyped(KeyEvent e) {
                e.consume();
            }

            @Override
            public void keyReleased(KeyEvent e) {
                e.consume();
            }

            @Override
            public void keyPressed(KeyEvent e) {
                e.consume();
            }
        });
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

    public void beginROI(JButton button) {
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

    // Set the focus on the command bar search box
    public void focusSearchBox() {
        contentPane.getCommandBar().focusSearchBox();
    }


    public void selectGenomeFromList(String genomeId) {
        contentPane.getCommandBar().selectGenome(genomeId);
    }


    public void defineGenome(javax.swing.ProgressMonitor monitor) {

        ProgressBar.ProgressDialog progressDialog = null;
        File archiveFile = null;

        try {
            GenomeBuilderDialog genomeBuilderDialog = new GenomeBuilderDialog(mainFrame, this);
            genomeBuilderDialog.setVisible(true);

            File genomeZipFile = genomeBuilderDialog.getArchiveFile();
            if (genomeBuilderDialog.isCanceled() || genomeZipFile == null) {
                return;
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
                monitor.setProgress(100);
            }

        } catch (MaximumContigGenomeException e) {

            String genomePath = "";
            if (archiveFile != null) {
                genomePath = archiveFile.getAbsolutePath();
            }

            log.error("Failed to define genome: " + genomePath, e);

            JOptionPane.showMessageDialog(mainFrame, "Failed to define genome " +
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
        }
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

        UIUtilities.invokeAndWaitOnEventThread(() -> {

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
        });
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


    final public void doViewPreferences() {

        // 2.x releases -- swing editor
        if (Globals.VERSION.contains("2.4")) {
            PreferencesEditor dialog = new PreferencesEditor(this.mainFrame, true);
            dialog.setVisible(true);
        } else {
            // 3.0 releases -- javafx
            try {
                PreferenceEditorFX.open(this.mainFrame);
            } catch (IOException e) {
                log.error("Error openining preference dialog", e);
            }
        }
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
            PreferencesManager.getPreferences().remove(RECENT_SESSIONS);
            PreferencesManager.getPreferences().setRecentSessions(recentSessions);
        }

    }

    final public void doShowAttributeDisplay(boolean enableAttributeView) {

        boolean oldState = PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY);

        // First store the newly requested state
        if (oldState != enableAttributeView) {
            PreferencesManager.getPreferences().setShowAttributeView(enableAttributeView);
            doRefresh();
        }


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
                    PreferencesManager.getPreferences().put(SHOW_DEFAULT_TRACK_ATTRIBUTES, true);
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

        if (extension == null) {
            extension = ".png";
            file = new File(file.getAbsolutePath() + extension);
        }

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

        File snapshotDirectory = PreferencesManager.getPreferences().getLastSnapshotDirectory();

        JFileChooser fc = new SnapshotFileChooser(snapshotDirectory, defaultFile);
        fc.showSaveDialog(mainFrame);
        File file = fc.getSelectedFile();

        // If a file selection was made
        if (file != null) {
            File directory = file.getParentFile();
            if (directory != null) {
                PreferencesManager.getPreferences().setLastSnapshotDirectory(directory);
            }

        }

        return file;
    }


    private void createZoomCursors() throws HeadlessException, IndexOutOfBoundsException {
        if (zoomInCursor == null || zoomOutCursor == null) {
            final Image zoomInImage = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_IN).getImage();
            final Image zoomOutImage = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_OUT).getImage();
            final Point hotspot = new Point(10, 10);
            zoomInCursor = createCustomCursor(zoomInImage, hotspot, "Zoom in", Cursor.CROSSHAIR_CURSOR);
            zoomOutCursor = createCustomCursor(zoomOutImage, hotspot, "Zoom out", Cursor.DEFAULT_CURSOR);

        }

    }


    private void createHandCursor() throws HeadlessException, IndexOutOfBoundsException {

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
                        fistCursor = createCustomCursor(handImage, new Point(8, 6), "Move", Cursor.HAND_CURSOR);
                        return false;
                    } else {
                        return true;
                    }
                }
            });
            if (ready) {
                try {
                    fistCursor = createCustomCursor(handImage, new Point(8, 6), "Move", Cursor.HAND_CURSOR);
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
                        dragNDropCursor = createCustomCursor(dragNDropImage, new Point(0, 0), "Drag and Drop", Cursor.CROSSHAIR_CURSOR);
                        return false;

                    } else {
                        return true;
                    }
                }
            });
            if (ready) {
                dragNDropCursor = createCustomCursor(dragNDropImage, new Point(0, 0), "Drag and Drop", Cursor.DEFAULT_CURSOR);
            }
        }
    }


    private Cursor createCustomCursor(Image image, Point hotspot, String name, int defaultCursor) {
        try {
            return mainFrame.getToolkit().createCustomCursor(image, hotspot, name);
        } catch (Exception e) {
            log.info("Could not create cursor: " + name);
            return Cursor.getPredefinedCursor(defaultCursor);
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

        List<Track> oldTracks = getAllTracks();

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

        contentPane.getMainPanel().resetPanels();

        //TODO -- this is a very blunt and dangerous way to clean up -- change to close files associated with this session
        SeekableFileStream.closeAllInstances();

        Set<Track> newTracks = new HashSet<>(getAllTracks());
        for (Track t : oldTracks) {

        }

        doRefresh();
        System.gc();
    }

    private void subscribeToEvents() {
        IGVEventBus.getInstance().subscribe(ViewChange.class, this);
        IGVEventBus.getInstance().subscribe(ShiftEvent.class, this);
        IGVEventBus.getInstance().subscribe(InsertionSelectionEvent.class, this);
        IGVEventBus.getInstance().subscribe(GenomeChangeEvent.class, this);
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

    public void setStatusBarMessag2(String message) {
        contentPane.getStatusBar().setMessage2(message);
    }

    public void setStatusBarMessage3(String message) {
        contentPane.getStatusBar().setMessage3(message);
    }

    public void enableStopButton(boolean enable) {
        contentPane.getStatusBar().enableStopButton(enable);
    }

    /**
     * Resets factory settings. this is not the same as reset user defaults
     * DO NOT DELETE used when debugging
     */
    public void resetToFactorySettings() {

        try {
            PreferencesManager.getPreferences().clear();
            boolean isShow = PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY);
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

        // check to see if any files in session file are on protected (oauth) server. If
        // so, make sure user is logged into
        // server before -proceeding
        OAuthUtils.checkServerLogin(sessionPath);

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
     * @see IGVPreferences#getPersistent(String, String)
     */
    public static String getPersistent(String key, String def) {
        if (IGV.hasInstance()) {
            return IGV.getInstance().getSession().getPersistent(key, def);
        } else {
            return PreferencesManager.getPreferences().getPersistent(key, def);
        }
    }


    /**
     * Reset the default status message, which is the number of tracks loaded.
     */
    public void resetStatusMessage() {
        contentPane.getStatusBar().setMessage("" +
                getVisibleTrackCount() + " tracks loaded");

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
        return contentPane != null && contentPane.getCommandBar().getDetailsBehavior() == ShowDetailsBehavior.CLICK;
    }

    public boolean isShowDetailsOnHover() {
        return contentPane != null && contentPane.getCommandBar().getDetailsBehavior() == ShowDetailsBehavior.HOVER;
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
     * <p>
     * NOTE: this must be run on the event tread as UI components are added here.
     *
     * @param locators
     */
    public void loadResources(Collection<ResourceLocator> locators) {

        log.info("Loading " + locators.size() + " resources.");
        final MessageCollection messages = new MessageCollection();

        for (final ResourceLocator locator : locators) {

            // If its a local file, check explicitly for existence (rather than rely on exception)
            if (locator.isLocal()) {
                File trackSetFile = new File(locator.getPath());
                if (!trackSetFile.exists()) {
                    messages.append("File not found: " + locator.getPath() + "\n");
                    continue;
                }
            }

            try {
                List<Track> tracks = load(locator);
                addTracks(tracks, locator);
            } catch (Exception e) {
                log.error("Error loading track", e);
                messages.append("Error loading " + locator + ": " + e.getMessage());
            }

        }
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

            }
        }

        return newTracks;
    }


    /**
     * Load the data file into the specified panel.   Triggered via drag and drop.
     */
    public void load(final ResourceLocator locator, final TrackPanel panel) throws DataLoadException {

        // If this is a session  TODO -- need better "is a session?" test
        if (locator.getPath().endsWith(".xml") || locator.getPath().endsWith(("session"))) {
            boolean merge = false;  // TODO -- ask user?
            this.doRestoreSession(locator.getPath(), null, merge);
        } else {
            // Not a session, load into target panel
            Runnable runnable = () -> {
                List<Track> tracks = load(locator);
                panel.addTracks(tracks);
                doRefresh();
            };
            LongRunningTask.submit(runnable);
        }
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
        } else if (PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
            return getTrackPanel(DATA_PANEL_NAME);
        } else if (TrackLoader.isAlignmentTrack(locator.getTypeString())) {
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


        if (isAnnotationFile(filename)) {
            return getTrackPanel(FEATURE_PANEL_NAME);
        } else {
            return getTrackPanel(DATA_PANEL_NAME);
        }
    }

    private boolean isAnnotationFile(String filename) {
        return filename.contains("refflat") || filename.contains("ucscgene") ||
                filename.contains("genepred") || filename.contains("ensgene") ||
                filename.contains("refgene") ||
                filename.endsWith("gff") || filename.endsWith("gtf") ||
                filename.endsWith("gff3") || filename.endsWith("embl") ||
                filename.endsWith("bed") || filename.endsWith("gistic") ||
                filename.endsWith("bedz") || filename.endsWith("repmask") ||
                filename.contains("dranger") || filename.endsWith("ucscsnp");
    }

    public void reset() {
        groupByAttribute = null;
        for (TrackPanel sp : getTrackPanels()) {
            if (DATA_PANEL_NAME.equals(sp.getName())) {
                sp.reset();
                break;
            }
        }

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
    public void groupAlignmentTracks(AlignmentTrack.GroupOption option, String tag, Range pos) {
        final IGVPreferences prefMgr = PreferencesManager.getPreferences();
        prefMgr.put(SAM_GROUP_OPTION, option.toString());
        if (option == AlignmentTrack.GroupOption.TAG && tag != null) {
            prefMgr.put(SAM_GROUP_BY_TAG, tag);
        }
        if (option == AlignmentTrack.GroupOption.BASE_AT_POS && pos != null) {
            prefMgr.put(SAM_GROUP_BY_POS, pos.getChr() + " " + String.valueOf(pos.getStart()));
        }
        for (Track t : getAllTracks()) {
            if (t instanceof AlignmentTrack) {
                ((AlignmentTrack) t).groupAlignments(option, tag, pos);
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
            Track t = allTracks.get(i);
            if (t.isVisible()) {
                t.setSelected(true);
            }
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
        removeTracks(tracksToRemove, true);
    }

    public void removeTracks(Collection<? extends Track> tracksToRemove, boolean dispose) {

        // Make copy of list as we will be modifying the original in the loop
        List<TrackPanel> panels = getTrackPanels();
        for (TrackPanel trackPanel : panels) {
            trackPanel.removeTracks(tracksToRemove);

            if (!trackPanel.hasTracks()) {
                removeDataPanel(trackPanel.getName());
            }
        }

        for (Track t : tracksToRemove) {
            if (t instanceof IGVEventObserver) {
                IGVEventBus.getInstance().unsubscribe((IGVEventObserver) t);
            }
        }

        if (dispose) {
            for (Track t : tracksToRemove) {
                t.dispose();
            }
        }
    }

    /**
     * Add gene and sequence tracks.  This is called upon switching genomes.
     *
     * @param newGeneTrack
     * @param
     */
    public void setGenomeTracks(FeatureTrack newGeneTrack) {

        TrackPanel panel = PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY) ?
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
        revalidateTrackPanels();
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
        IGVEventBus.getInstance().post(new TrackGroupEvent(this));
    }


    private void resetGroups() {
        log.debug("Resetting Groups");
        for (TrackPanel trackPanel : getTrackPanels()) {
            trackPanel.groupTracksByAttribute(groupByAttribute);
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
        ProgressBar.ProgressDialog progressDialog;

        StartupRunnable(Main.IGVArgs args) {
            this.igvArgs = args;

        }

        @Override
        public void run() {

            final boolean runningBatch = igvArgs.getBatchFile() != null;
            BatchRunner.setIsBatchMode(runningBatch);


            UIUtilities.invokeOnEventThread(() -> mainFrame.setIconImage(getIconImage()));
            if (Globals.IS_MAC) {
                setAppleDockIcon();
            }

            final IGVPreferences preferenceManager = PreferencesManager.getPreferences();

            boolean genomeLoaded = false;
            if (igvArgs.getGenomeId() != null) {
                String genomeId = igvArgs.getGenomeId();
                try {
                    GenomeManager.getInstance().loadGenomeById(genomeId);
                    genomeLoaded = true;
                } catch (IOException e) {
                    MessageUtils.showErrorMessage("Error loading genome: " + genomeId, e);
                    log.error("Error loading genome: " + genomeId, e);
                }
            }


            if (genomeLoaded == false && igvArgs.getSessionFile() == null) {
                String genomeId = preferenceManager.getDefaultGenome();
                try {
                    GenomeManager.getInstance().loadGenomeById(genomeId);
                    genomeLoaded = true;
                } catch (Exception e) {
                    MessageUtils.showErrorMessage("Error loading genome: " + genomeId, e);
                    log.error("Error loading genome: " + genomeId, e);
                }
            }

            if (genomeLoaded == false && igvArgs.getSessionFile() == null) {
                String genomeId = GenomeManager.DEFAULT_GENOME.getId();
                try {
                    GenomeManager.getInstance().loadGenomeById(genomeId);
                } catch (IOException e) {
                    MessageUtils.showErrorMessage("Error loading genome: " + genomeId, e);
                    log.error("Error loading genome: " + genomeId, e);
                }
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
                    String decodedString = igvArgs.getDataFileString();
                    String[] dataFiles = decodedString.split(",");
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
            boolean portEnabled = preferenceManager.getAsBoolean(PORT_ENABLED);
            String portString = igvArgs.getPort();
            if (portEnabled || portString != null) {
                // Command listener thread
                int port = preferenceManager.getAsInt(PORT_NUMBER);
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
                        Globals.setBatch(false);   // Set to false for startup execution -- otherwise we block the event thread
                        LongRunningTask.submit(new BatchRunner(igvArgs.getBatchFile()));
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
                DesktopIntegration.setDockIcon(image);
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
            pluginClassNames.addAll(Arrays.asList(PreferencesManager.getPreferences().getIGVPluginList()));
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

    public static void copySequenceToClipboard(Genome genome, String chr, int start, int end, Strand strand) {
        try {
            IGV.getMainFrame().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
            byte[] seqBytes = genome.getSequence(chr, start, end);

            if (seqBytes == null) {
                MessageUtils.showMessage("Sequence not available");
            } else {
                String sequence = new String(seqBytes);

                SequenceTrack sequenceTrack = IGV.getInstance().getSequenceTrack();
                if (strand == Strand.NEGATIVE || (sequenceTrack != null && sequenceTrack.getStrand() == Strand.NEGATIVE)) {
                    sequence = SequenceTrack.getReverseComplement(sequence);
                }
                StringUtils.copyTextToClipboard(sequence);
            }

        } finally {
            IGV.getMainFrame().setCursor(Cursor.getDefaultCursor());
        }
    }


    /**
     * Wrapper for igv.wait(timeout).   Used during unit tests.
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

    public void receiveEvent(Object event) {

        if (event instanceof ViewChange || event instanceof InsertionSelectionEvent) {
            revalidateTrackPanels();   // TODO -- this seems extreme
        } else if (event instanceof ShiftEvent) {
            revalidateTrackPanels();
        } else if (event instanceof GenomeChangeEvent) {
            doRefresh();
        } else {
            log.info("Unknown event type: " + event.getClass());
        }

    }


    public void repaint() {
        mainFrame.repaint();
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

    final public void doRefresh() {
        contentPane.getMainPanel().revalidate();
        mainFrame.repaint();
        getContentPane().repaint();
    }

    /**
     * Repaint the header and data panels.
     * <p/>
     * Note:  If running in Batch mode we force synchronous painting.  This is necessary as the
     * paint() command triggers loading of data.  If allowed to proceed asynchronously the "snapshot" batch command
     * might execute before the data from a previous command has loaded.
     */
    public void revalidateTrackPanels() {

        UIUtilities.invokeOnEventThread(() -> {

            if (Globals.isBatch()) {
                contentPane.revalidateTrackPanels();
                rootPane.paintImmediately(rootPane.getBounds());
            } else {
                contentPane.revalidateTrackPanels();
                rootPane.repaint();

            }
        });
    }

    public void repaintNamePanels() {
        for (TrackPanel tp : getTrackPanels()) {
            tp.getScrollPane().getNamePanel().repaint();
        }
    }

//
//    NOTE:  MAC ONLY,  WILL NOT COMPILE ON OTHER PLATFORMS
//    private void getRealDPI() {
//        // find the display device of interest
//        final GraphicsDevice defaultScreenDevice = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice();
//
//        // on OS X, it would be CGraphicsDevice
//        if (defaultScreenDevice instanceof CGraphicsDevice) {
//            final CGraphicsDevice device = (CGraphicsDevice) defaultScreenDevice;
//
//            // this is the missing correction factor, it's equal to 2 on HiDPI a.k.a. Retina displays
//            final int scaleFactor = device.getScaleFactor();
//
//            // now we can compute the real DPI of the screen
//            final double realDPI = scaleFactor * (device.getXResolution() + device.getYResolution()) / 2;
//
//            System.out.println("scale factor = " + scaleFactor + "    realDPI = " + realDPI);
//        }
//    }


}
