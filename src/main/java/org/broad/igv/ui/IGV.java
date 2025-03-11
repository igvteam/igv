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
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.annotations.ForTesting;
import org.broad.igv.batch.BatchRunner;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.event.*;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.*;
import org.broad.igv.lists.GeneList;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesEditor;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.InsertionSelectionEvent;
import org.broad.igv.sam.SortOption;
import org.broad.igv.session.*;
import org.broad.igv.session.autosave.AutosaveTimerTask;
import org.broad.igv.session.autosave.SessionAutosaveManager;
import org.broad.igv.track.*;
import org.broad.igv.ui.WaitCursorManager.CursorToken;
import org.broad.igv.ui.dnd.GhostGlassPane;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.*;
import org.broad.igv.util.*;
import org.broad.igv.variant.VariantTrack;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.List;
import java.util.*;
import java.util.Timer;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.prefs.Preferences;
import java.util.stream.Collectors;

import static org.broad.igv.prefs.Constants.*;

/**
 * Represents an IGV instance, consisting of a main window and associated model.
 *
 * @author jrobinso
 */
public class IGV implements IGVEventObserver {

    private static Logger log = LogManager.getLogger(IGV.class);
    private static IGV theInstance;

    public static final String DATA_PANEL_NAME = "DataPanel";
    public static final String FEATURE_PANEL_NAME = "FeaturePanel";

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

    /**
     * Object to hold state that defines a user session.  There is always a user session, even if not initialized
     * from a "session" file.
     */
    private Session session;

    /**
     * Timer for triggering periodic autosave of current session
     */
    private Timer sessionAutosaveTimer = new Timer();

    // Misc state
    private Map<String, List<Track>> overlayTracksMap = new HashMap<>();
    private Set<Track> overlaidTracks = new HashSet<>();
    private RecentFileSet recentSessionList;
    private RecentUrlsSet recentUrlsList;

    // Vertical line that follows the mouse
    private boolean rulerEnabled;

    public static IGV createInstance(Frame frame, Main.IGVArgs igvArgs) {
        if (theInstance != null) {
            throw new RuntimeException("Only a single instance is allowed.");
        }
        theInstance = new IGV(frame, igvArgs);
        return theInstance;
    }

    public static IGV getInstance() {
        if (theInstance == null) {
            throw new RuntimeException("IGV has not been initialized.  Must call createInstance(Frame) first");
        }
        return theInstance;
    }

    public static boolean hasInstance() {
        return theInstance != null;
    }

    @ForTesting
    static void destroyInstance() {
        IGVMenuBar.destroyInstance();
        theInstance = null;
    }

    /**
     * Creates new IGV
     */
    private IGV(Frame frame, Main.IGVArgs igvArgs) {

        theInstance = this;

        final IGVPreferences preferences = PreferencesManager.getPreferences();

        session = new Session(null);

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
        // get the current main screen bounds
        Dimension screenBounds = Toolkit.getDefaultToolkit().getScreenSize();
        Rectangle applicationBounds = preferences.getApplicationFrameBounds();

        // get info of all screens and store them into the array
        GraphicsEnvironment graphEnv = GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice[] graphDev = graphEnv.getScreenDevices();
        Rectangle[] boundsArr = new Rectangle[graphDev.length];

        for (int i = 0; i < graphDev.length; ++i) {
            GraphicsConfiguration curCon = graphDev[i].getDefaultConfiguration();
            boundsArr[i] = curCon.getBounds();
        }

        //set a flag which indicates if the user preference is empty, or if the (x,y) in the user preference is not contained in any screen
        //default is empty or not contained
        boolean isNullOrNotContained = true;

        if (applicationBounds != null) {
            //Iterate over each screen value to find if there is currently a screen that can contain these values.
            int userX = applicationBounds.x;
            int userY = applicationBounds.y;
            double userMaxX = applicationBounds.getMaxX();
            double userMaxY = applicationBounds.getMaxY();
            for (Rectangle curScreen : boundsArr) {
                if (curScreen.contains(userX, userY)) {
                    isNullOrNotContained = false;
                    if (userMaxX >= curScreen.getMaxX() || userMaxY >= curScreen.getMaxY()) {
                        applicationBounds = new Rectangle(curScreen.x, curScreen.y, Math.min(1150, curScreen.width), Math.min(800, curScreen.height));
                    }
                    break;
                }
            }
        }
        if (isNullOrNotContained) {
            // user's preference is null or the (x,y) in user's preference is not contained in any screen
            // set the application to the main screen
            applicationBounds = new Rectangle(0, 0, Math.min(1150, screenBounds.width), Math.min(800, screenBounds.height));
        }
        mainFrame.setBounds(applicationBounds);

        subscribeToEvents();


        // Start running periodic autosaves (unless the user has specified not to retain timed autosaves)

        if (PreferencesManager.getPreferences().getAsInt(Constants.AUTOSAVES_TO_KEEP) > 0) {
            int timerDelay = PreferencesManager.getPreferences().getAsInt(AUTOSAVE_FREQUENCY) * 60000; // Convert timer delay to ms
            sessionAutosaveTimer.scheduleAtFixedRate(new AutosaveTimerTask(this), timerDelay, timerDelay);
        }
    }


    public JRootPane getRootPane() {
        return rootPane;
    }

    public Frame getMainFrame() {
        return mainFrame;
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
        repaint();
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

        Future toRet = null;
        if (locators != null && !locators.isEmpty()) {

            contentPane.getStatusBar().setMessage("Loading ...");

            NamedRunnable runnable = new NamedRunnable() {

                public void run() {
                    //Collect size statistics before loading
                    List<Map<TrackPanelScrollPane, Integer>> trackPanelAttrs = getTrackPanelAttrs();

                    final MessageCollection messages = new MessageCollection();
                    for (final ResourceLocator locator : locators) {

                        // If it's a local file, check explicitly for existence (rather than rely on exception)
                        if (locator.isLocal()) {
                            File trackSetFile = new File(locator.getPath());
                            if (!trackSetFile.exists()) {
                                messages.append("File not found: " + locator.getPath() + "\n");
                                continue;
                            }
                        }

                        try {
                            List<Track> tracks = load(locator);
                            addTracks(tracks);

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


                    resetPanelHeights(trackPanelAttrs.get(0), trackPanelAttrs.get(1));
                    showLoadedTrackCount();
                    revalidateTrackPanels();
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
        try {
            PreferencesEditor.open(this.mainFrame);
        } catch (Exception e) {
            log.error("Error opening preference dialog", e);
        }
    }

    final public void saveStateForExit() {

        // Store recent sessions
        RecentFileSet recentSessions = getRecentSessionList();
        if (!recentSessions.isEmpty()) {
            PreferencesManager.getPreferences().setRecentSessions(recentSessions.asString());
        }

        // Store recent files
        RecentUrlsSet recentUrls = getRecentUrls();
        if (!recentUrls.isEmpty()) {
            PreferencesManager.getPreferences().setRecentUrls(recentUrls.asString());
        }

        // Stop the timer that is triggering the timed autosave
        this.stopTimedAutosave();

        // Autosave current session if configured to do so
        if (PreferencesManager.getPreferences().getAsBoolean(AUTOSAVE_ON_EXIT)) {
            try {
                SessionAutosaveManager.saveExitSessionAutosaveFile(session);
            } catch (Exception e) {
                log.error("Error autosaving session", e);
            }
        }


    }

    final public void doShowAttributeDisplay(boolean enableAttributeView) {

        boolean oldState = PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY);

        // First store the newly requested state
        if (oldState != enableAttributeView) {
            PreferencesManager.getPreferences().setShowAttributeView(enableAttributeView);
            repaint();
        }
    }


    /**
     * Set the attributes to show in the attribute panel for this session.
     */
    final public void doSelectDisplayableAttribute() {

        List<String> allAttributes = AttributeManager.getInstance().getAttributeNames();
        Set<String> hiddenAttributes = IGV.getInstance().getSession().getHiddenAttributes();
        final CheckListDialog dlg = new CheckListDialog(mainFrame, allAttributes, hiddenAttributes, false);
        dlg.setVisible(true);

        if (!dlg.isCanceled()) {
            IGV.getInstance().getSession().setHiddenAttributes(dlg.getNonSelections());
            revalidateTrackPanels();
        }
    }


    final public void saveImage(Component target, String extension) {
        saveImage(target, "igv_snapshot", extension);
    }

    final public void saveImage(Component target, String title, String extension) {
        if ("png".equalsIgnoreCase(extension) || "svg".equalsIgnoreCase(extension)) {
            contentPane.getStatusBar().setMessage("Creating image...");
            File defaultFile = new File(title + "." + extension);
            createSnapshot(target, defaultFile);
        }
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
        } catch (Exception e) {
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
     * @throws IOException
     * @api
     * @see ImageFileTypes.Type
     */
    public String createSnapshotNonInteractive(Component target, File file, boolean batch) throws Exception {

        log.debug("Creating snapshot: " + file.getName());

        String extension = FileUtils.getFileExtension(file.getAbsolutePath());

        if (extension == null) {
            extension = ".png";
            file = new File(file.getAbsolutePath() + extension);
        }

        ImageFileTypes.Type type = ImageFileTypes.getImageFileType(extension);

        if (type == ImageFileTypes.Type.NULL) {
            String message = "ERROR: Unknown file extension " + extension;
            log.error(message);
            return message;
        } else if (type == ImageFileTypes.Type.EPS || type == ImageFileTypes.Type.JPEG) {
            String message = "ERROR: " + type + " output is not supported.  Try '.png' or '.svg'";
            log.error(message);
            return message;
        }

        try {
            return SnapshotUtilities.doComponentSnapshot(target, file, type, batch);
        } finally {
            log.debug("Finished creating snapshot: " + file.getName());
        }
    }

    public File selectSnapshotFile(File defaultFile) {

        File snapshotDirectory = PreferencesManager.getPreferences().getLastSnapshotDirectory();

        // JFileChooser fc = new SnapshotFileChooser(snapshotDirectory, defaultFile);
        FileDialog fc = new FileDialog(mainFrame, "Save image", FileDialog.SAVE);
        if (snapshotDirectory != null) {
            fc.setDirectory(snapshotDirectory.getAbsolutePath());
        }
        fc.setFile(defaultFile.getName());
        fc.setFilenameFilter((dir, name) ->
                name.endsWith(".jpeg") || name.endsWith(".jpg") || name.endsWith(".png") || name.endsWith(".svg"));
        fc.setVisible(true);
        String file = fc.getFile();
        // If a file selection was made
        if (file != null) {
            String directory = fc.getDirectory();
            if (directory != null) {
                PreferencesManager.getPreferences().setLastSnapshotDirectory(directory);
            }
            return new File(directory, file);
        } else {
            return null;
        }
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
                    log.warn("Warning: could not create fistCursor");
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
            log.warn("Could not create cursor: " + name);
            return Cursor.getPredefinedCursor(defaultCursor);
        }
    }

    private void subscribeToEvents() {
        IGVEventBus.getInstance().subscribe(ViewChange.class, this);
        IGVEventBus.getInstance().subscribe(InsertionSelectionEvent.class, this);
        IGVEventBus.getInstance().subscribe(GenomeChangeEvent.class, this);
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
            repaint();

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


    /**
     * Scroll panel(s) containing track with specified name to the top.  Supports batch command "scrolltotrack" *
     *
     * @param trackName
     * @return
     */
    public boolean scrollToTrack(String trackName) {
        boolean found = false;
        for (TrackPanel tp : getTrackPanels()) {
            if (tp.getScrollPane().getNamePanel().scrollTo(trackName)) {
                found = true;
            }
        }
        return found;
    }


    /**
     * Scroll all panels to the top (position 0).  Supports batch command "scrolltotop"
     *
     * @return
     */
    public void scrollToTop() {
        for (TrackPanel tp : getTrackPanels()) {
            tp.getScrollPane().getNamePanel().scrollToPosition(0);

        }
    }


    public Session getSession() {
        return session;
    }


    /**
     * Reset session state, and associate the session object with the given path to a session file.
     *
     * @param sessionPath
     */
    public void resetSession(String sessionPath) {

        session.reset(sessionPath);
        AttributeManager.getInstance().clearAllAttributes();
        mainFrame.setTitle(sessionPath == null ? UIConstants.APPLICATION_NAME : sessionPath);
        menuBar.resetSessionActions();

        getMainPanel().resetPanels();   // Also clears all tracks
        getMainPanel().updatePanelDimensions();
        revalidateTrackPanels();
    }


    /**
     * Creates a new IGV session
     */
    public void newSession() {
        resetSession(null);
        Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
        if (currentGenome != null) {
            GenomeManager.getInstance().restoreGenomeTracks(currentGenome);
        }
        this.menuBar.disableReloadSession();
        goToLocus(GenomeManager.getInstance().getCurrentGenome().getHomeChromosome());
        revalidateTrackPanels();
    }


    /**
     * Load a session file, then jump to the specified locus if supplied.  This runs in the current thread, and should
     * not be called from the event dispatch thread.
     *
     * @param sessionPath
     * @param locus
     * @return true if successful
     */
    public boolean loadSession(String sessionPath, String locus) {

        InputStream inputStream = null;
        try {
            setStatusBarMessage("Opening session...");

            inputStream = new BufferedInputStream(ParsingUtils.openInputStreamGZ(new ResourceLocator(sessionPath)));
            boolean success = loadSessionFromStream(sessionPath, inputStream);

            if (success) {

                session.setPath(sessionPath);

                String searchText = locus == null ? session.getLocus() : locus;
                if (!FrameManager.isGeneListMode() && searchText != null && searchText.trim().length() > 0) {
                    goToLocus(searchText);
                }

                mainFrame.setTitle(UIConstants.APPLICATION_NAME + " - Session: " + sessionPath);

                getRecentSessionList().add(sessionPath);
                this.menuBar.enableReloadSession();

                //If there's a RegionNavigatorDialog, kill it.
                //this could be done through the Observer that RND uses, I suppose.  Not sure that's cleaner
                RegionNavigatorDialog.destroyInstance();
            }

            return success;

        } catch (Exception e) {
            String message = "Error loading session session: " + e.getMessage();
            MessageUtils.showMessage(message);
            getRecentSessionList().remove(sessionPath);
            log.error(e);
            return false;
        } finally {
            if (inputStream != null) {
                try {
                    inputStream.close();
                } catch (IOException iOException) {
                    log.error("Error closing session stream", iOException);
                }
            }
            resetStatusMessage();
        }
    }

    /**
     * Load a session from the input stream.  Currently there are
     * 2 sources for the input stream (1) a local or remote file, and (2) an in memory byte array.  The latter is
     * used to support the "reloadTracks" menu function.
     *
     * @param sessionPath
     * @param inputStream
     * @return
     * @throws IOException
     */

    public boolean loadSessionFromStream(String sessionPath, InputStream inputStream) throws IOException {

        final SessionReader sessionReader;
        if (sessionPath != null && (sessionPath.endsWith(".session") || sessionPath.endsWith(".session.txt"))) {
            sessionReader = new UCSCSessionReader(this);
        } else if (sessionPath != null && (sessionPath.endsWith(".idxsession") || sessionPath.endsWith(".idxsession.txt"))) {
            sessionReader = new IndexAwareSessionReader(this);
        } else {
            sessionReader = new IGVSessionReader(this);
        }

        sessionReader.loadSession(inputStream, session, sessionPath);

        double[] dividerFractions = session.getDividerFractions();
        if (dividerFractions != null) {
            contentPane.getMainPanel().setDividerFractions(dividerFractions);
        }
        session.clearDividerLocations();

        revalidateTrackPanels();
        return true;
    }

    /**
     * Saves current session to {@code targetFile}. As a side effect,
     * sets the current sessions path (does NOT set the last session directory)
     *
     * @param targetFile
     * @throws IOException
     */
    public void saveSession(File targetFile) throws IOException {
        (new SessionWriter()).saveSession(session, targetFile);

        String sessionPath = targetFile.getAbsolutePath();
        session.setPath(sessionPath);
        mainFrame.setTitle(UIConstants.APPLICATION_NAME + " - Session: " + sessionPath);

        getRecentSessionList().add(sessionPath);
        this.menuBar.enableReloadSession();

        // No errors so save last location
        PreferencesManager.getPreferences().setLastTrackDirectory(targetFile.getParentFile());

    }

    /**
     * Reset the default status message, which is the number of tracks loaded.
     */
    public void resetStatusMessage() {
        UIUtilities.invokeAndWaitOnEventThread(() -> contentPane.getStatusBar().setMessage("" +
                getVisibleTrackCount() + " tracks loaded"));

    }

    public void showLoadedTrackCount() {

        final int visibleTrackCount = getVisibleTrackCount();
        contentPane.getStatusBar().setMessage("" +
                visibleTrackCount + (visibleTrackCount == 1 ? " track" : " tracks"));
    }

    private void closeWindow(final ProgressBar.ProgressDialog progressDialog) {
        UIUtilities.invokeOnEventThread(() -> progressDialog.setVisible(false));
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

    public void tweakPanelDivider() {
        contentPane.getMainPanel().tweakPanelDivider();
    }

    public void removeDataPanel(String name) {
        contentPane.getMainPanel().removeDataPanel(name);
    }

    public void removeTrackPanel(TrackPanel trackPanel) {
        contentPane.getMainPanel().removeTrackPanel(trackPanel);
    }

    public boolean panelIsRemovable(TrackPanel trackPanel) {
        return contentPane.getMainPanel().panelIsRemovable(trackPanel);
    }

    public MainPanel getMainPanel() {
        return contentPane.getMainPanel();
    }

    public RecentFileSet getRecentSessionList() {
        if (recentSessionList == null) {
            recentSessionList = PreferencesManager.getPreferences().getRecentSessions();
            //remove sessions that no longer exist
            recentSessionList.removeIf(file -> !(new File(file)).exists());
        }
        return recentSessionList;
    }

    public RecentUrlsSet getRecentUrls() {
        if (recentUrlsList == null) {
            recentUrlsList = PreferencesManager.getPreferences().getRecentUrls();
        }
        return recentUrlsList;
    }

    /**
     * Add new values to the recent URLS set.  Calling this method rather than adding them directly
     * allows showing the menu when the first URL is added to the collection.
     * @param toAdd
     */
    public void addToRecentUrls(Collection<ResourceLocator> toAdd) {
        RecentUrlsSet recentFiles = getRecentUrls();
        recentFiles.addAll(toAdd);
        if (!recentFiles.isEmpty()) {
            menuBar.showRecentFilesMenu();
        }
    }


    public IGVContentPane getContentPane() {
        return contentPane;
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
     * Add the specified tracks to the appropriate panel. The list of tracks share a common file (locator).  Panel
     * is chosen based on characteristics of the {@code locator}.
     */
    public void addTracks(List<Track> trackList) {

        // Group tracks by locator.
        Map<String, List<Track>> map = trackList.stream().collect(Collectors.groupingBy(t -> t.getResourceLocator().getPath()));

        for(List<Track> tracks : map.values()) {
            Track representativeTrack = tracks.get(0);
            TrackPanel panel = getPanelFor(representativeTrack);
            panel.addTracks(tracks);
        }
    }


    /**
     * Load a resource and return the tracks.
     * Does not add tracks to igv instance
     *
     * @param locator
     * @return A list of loaded tracks
     */
    public List<Track> load(ResourceLocator locator) throws DataLoadException {

        try {
            if (IGV.hasInstance()) {
                IGV.getInstance().setStatusBarMessage3("Loading " + locator.getPath());
            }

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

                    TrackProperties properties = locator.getTrackProperties();
                    if (properties != null) {
                        track.setProperties(properties);
                    }
                }
            }
            //revalidateTrackPanels();
            return newTracks;
        } finally {
            if (IGV.hasInstance()) {
                IGV.getInstance().setStatusBarMessage3("");
            }
        }
    }


    /**
     * Load the data file into the specified panel.   Triggered via drag and drop.
     */
    public void load(final ResourceLocator locator, final TrackPanel panel) throws DataLoadException {

        // If this is a session  TODO -- need better "is a session?" test
        if (SessionReader.isSessionFile(locator.getPath())) {
            LongRunningTask.submit(() -> this.loadSession(locator.getPath(), null));
        } else {
            // Not a session, load into target panel
            Runnable runnable = () -> {
                List<Track> tracks = load(locator);
                panel.addTracks(tracks);
                repaint();
            };
            LongRunningTask.submit(runnable);
        }
    }

    /**
     * Return a DataPanel for the given track.
     *
     * @param track
     * @return
     */
    public TrackPanel getPanelFor(Track track) {

        if (PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
            return getTrackPanel(DATA_PANEL_NAME);
        }

        ResourceLocator locator = track.getResourceLocator();
        final String format = locator.getFormat();
        if (locator.getPanelName() != null) {
            return getTrackPanel(locator.getPanelName());
        } else if ("alist".equals(format)) {
            return getVcfBamPanel();
        } else if (PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
            return getTrackPanel(DATA_PANEL_NAME);
        } else if (TrackLoader.isAlignmentTrack(format)) {
            String newPanelName = "Panel" + System.currentTimeMillis();
            return addDataPanel(newPanelName).getTrackPanel();
        } else  if (track instanceof VariantTrack && ((VariantTrack) track).getAllSamples().size() > 10) {
            String newPanelName = "Panel" + System.currentTimeMillis();
            return addDataPanel(newPanelName).getTrackPanel();
        } else {
            // Default
            return getTrackPanel(DATA_PANEL_NAME);
        }
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

    private boolean isAnnotationFile(String format) {
        Set<String> annotationFormats = new HashSet<>(Arrays.asList("refflat", "ucscgene",
                "genepred", "ensgene", "refgene", "gff", "gtf", "gff3", "embl", "bed", "gistic",
                "bedz", "repmask", "dranger", "ucscsnp", "genepredext", "bigbed"));
        return annotationFormats.contains(format);
    }


    public void sortAlignmentTracks(SortOption option, String tag, final boolean invertSort) {
        sortAlignmentTracks(option, null, tag, invertSort, null);
    }

    public void sortAlignmentTracks(SortOption option, Double location, String tag, boolean invertSort, Set<String> priorityRecords) {
        List<AlignmentTrack> alignmentTracks = getAllTracks().stream()
                .filter(track -> track instanceof AlignmentTrack)
                .map(track -> (AlignmentTrack) track)
                .peek(track -> track.sortRows(option, location, tag, invertSort, priorityRecords))
                .collect(Collectors.toList());
        this.repaint(alignmentTracks);
    }

    /**
     * Group all alignment tracks by the specified option.
     *
     * @param option
     * @api
     */
    public void groupAlignmentTracks(AlignmentTrack.GroupOption option, String tag, Range pos) {

        List<Track> alignmentTracks = getAllTracks().stream()
                .filter(track -> track instanceof AlignmentTrack)
                .collect(Collectors.toList());
        for (Track t : alignmentTracks) {
            ((AlignmentTrack) t).groupAlignments(option, tag, pos);
        }
        this.repaint(alignmentTracks);
    }

    /**
     * Group all alignment tracks by the specified option.
     *
     * @param option
     * @api
     */
    public void colorAlignmentTracks(AlignmentTrack.ColorOption option, String tag) {

        List<Track> alignmentTracks = getAllTracks().stream()
                .filter(track -> track instanceof AlignmentTrack)
                .collect(Collectors.toList());
        for (Track t : alignmentTracks) {
            final AlignmentTrack alignmentTrack = (AlignmentTrack) t;
            alignmentTrack.setColorOption(option);
            if (option == AlignmentTrack.ColorOption.BISULFITE && tag != null) {
                try {
                    AlignmentTrack.BisulfiteContext context = AlignmentTrack.BisulfiteContext.valueOf(tag);
                    alignmentTrack.setBisulfiteContext(context);
                } catch (IllegalArgumentException e) {
                    log.error("Error setting bisulfite context for: " + tag, e);
                }
            } else if (tag != null) {
                alignmentTrack.setColorByTag(tag);
            }
        }
        this.repaint(alignmentTracks);
    }

    public void packAlignmentTracks() {
        for (Track t : getAllTracks()) {
            if (t instanceof AlignmentTrack) {
                ((AlignmentTrack) t).packAlignments();
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
        return Lists.newArrayList(Iterables.filter(getAllTracks(), FeatureTrack.class));
    }

    public List<DataTrack> getDataTracks() {
        return Lists.newArrayList(Iterables.filter(getAllTracks(), DataTrack.class));
    }

    public List<AlignmentTrack> getAlignmentTracks() {
        return Lists.newArrayList(Iterables.filter(getAllTracks(), AlignmentTrack.class));
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


    public Set<TrackType> getLoadedTypes() {
        Set<TrackType> types = new HashSet<>();
        for (Track t : getAllTracks()) {
            if (t != null) {
                TrackType type = t.getTrackType();
                types.add(type);
            }
        }
        return types;
    }


    public Set<String> getLoadedPaths() {
        return getDataResourceLocators().stream()
                .map(rl -> rl.getPath())
                .collect(Collectors.toSet());
    }


    public void setAllTrackHeights(int newHeight) {
        for (Track track : getAllTracks()) {
            track.setHeight(newHeight, true);
        }

    }

    /**
     * Remove and dispose of tracks.  Removed tracks will not be usable afterwards.
     *
     * @param tracksToRemove
     */
    public void deleteTracks(Collection<? extends Track> tracksToRemove) {

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
            t.unload();
        }
        revalidateTrackPanels();
    }

    /**
     * Add gene and sequence tracks.  This is called upon switching genomes.
     *
     * @param
     */
    public void setSequenceTrack() {

        TrackPanel panel = PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY) ?
                getTrackPanel(DATA_PANEL_NAME) : getTrackPanel(FEATURE_PANEL_NAME);
        SequenceTrack newSeqTrack = new SequenceTrack("Reference sequence");
        panel.addTrack(newSeqTrack);

//        if (newGeneTrack != null) {
//            newGeneTrack.setAttributeValue(Globals.TRACK_NAME_ATTRIBUTE, newGeneTrack.getName());
//            newGeneTrack.setAttributeValue(Globals.TRACK_DATA_FILE_ATTRIBUTE, "");
//            newGeneTrack.setAttributeValue(Globals.TRACK_DATA_TYPE_ATTRIBUTE, newGeneTrack.getTrackType().toString());
//            panel.addTrack(newGeneTrack);
//        }
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
        repaint();
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

            Comparator<Track> c = (t1, t2) -> {
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

            };
            Collections.sort(tracks, c);

        }
    }


    ////////////////////////////////////////////////////////////////////////////////////////
    // Groups
    public String getGroupByAttribute() {
        return session.getGroupByAttribute();
    }


    public void setGroupByAttribute(String attributeName) {
        session.setGroupByAttribute(attributeName);
        resetGroups();
        // Some tracks need to respond to changes in grouping, fire notification event
        IGVEventBus.getInstance().post(new TrackGroupEvent());
    }


    private void resetGroups() {
        log.debug("Resetting Groups");
        for (TrackPanel trackPanel : getTrackPanels()) {
            trackPanel.groupTracksByAttribute(session.getGroupByAttribute());
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
     * Enables command port early, otherwise private URLs pointing to custom genomes cannot be accessed.
     * This is because CommandListener (http://localhost:65301) is needed for OAuth's redirect parameter.
     *
     * @param igvArgs: Used to specify a different port.
     */
    private static void startCommandsServer(Main.IGVArgs igvArgs, IGVPreferences prefMgr) {
        // Port # can be overriden with "-p" command line switch
        boolean portEnabled = prefMgr.getAsBoolean(PORT_ENABLED);
        String portString = igvArgs.getPort();

        if (portEnabled || portString != null) {
            // Command listener thread
            int port = prefMgr.getAsInt(PORT_NUMBER);
            if (portString != null) {
                port = Integer.parseInt(portString);
            }
            CommandListener.start(port);
        }
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

            final IGVPreferences preferences = PreferencesManager.getPreferences();

            // Start CommandsServer **before** loading the initial genome, as credentials might need to be set for
            // privately hosted genomes.
            startCommandsServer(igvArgs, preferences);

            UIUtilities.invokeAndWaitOnEventThread(() -> {
                mainFrame.setIconImage(getIconImage());
                if (Globals.IS_MAC) {
                    setAppleDockIcon();
                }
                mainFrame.setVisible(true);
            });

            // Load the initial genome.
            final boolean runningBatch = igvArgs.getBatchFile() != null;

            if (runningBatch) {

                BatchRunner.setIsBatchMode(true);
                try {
                    UIUtilities.invokeAndWaitOnEventThread(() -> {
                        String genomeId = preferences.getDefaultGenome();
                        BatchRunner batchRunner = (new BatchRunner(igvArgs.getBatchFile(), IGV.this));
                        batchRunner.runWithDefaultGenome(genomeId);
                    });
                } finally {
                    BatchRunner.setIsBatchMode(false);
                }

            } else {
                // Check whether autosave is set to load and exists
                boolean autosavePresent = false;
                try {
                    autosavePresent = SessionAutosaveManager.getMostRecentAutosaveFile().isPresent();
                } catch (Exception e) {
                    log.error("Failure trying to get most recent autosave file", e);
                }
                boolean loadAutosave = autosavePresent && PreferencesManager.getPreferences().getAsBoolean(AUTOLOAD_LAST_AUTOSAVE);

                boolean genomeLoaded = false;
                if (igvArgs.getGenomeId() != null) {
                    String genomeId = igvArgs.getGenomeId();
                    try {
                        genomeLoaded = GenomeManager.getInstance().loadGenomeById(genomeId);
                    } catch (IOException e) {
                        MessageUtils.showErrorMessage("Error loading genome: " + genomeId, e);
                        log.error("Error loading genome: " + genomeId, e);
                    }
                }
                // If we're not loading a session file, attempt to load a default genome file
                if (igvArgs.getSessionFile() == null && !loadAutosave && !genomeLoaded) {
                    String genomeId = preferences.getDefaultGenome();
                    try {
                        genomeLoaded = GenomeManager.getInstance().loadGenomeById(genomeId);
                    } catch (Exception e) {
                        MessageUtils.showErrorMessage("Error loading genome " + genomeId + "<br/>" + e.getMessage(), e);
                        genomeLoaded = false;
                    }

                    if (!genomeLoaded) {
                        Genome genome = Genome.nullGenome();
                        GenomeManager.getInstance().setCurrentGenome(genome);
                    }
                }

                if (igvArgs.getSessionFile() != null || igvArgs.getDataFileStrings() != null || loadAutosave) {

                    if (log.isDebugEnabled()) {
                        log.debug("Loading session data");
                    }

                    if (log.isDebugEnabled()) {
                        log.debug("Calling restore session");
                    }


                    if (igvArgs.getSessionFile() != null) {
                        IGV.getInstance().setStatusBarMessage3("Loading " + igvArgs.getSessionFile());
                        boolean success;
                        try {
                            success = false;
                            if (HttpUtils.isRemoteURL(igvArgs.getSessionFile())) {
                                boolean merge = false;
                                success = loadSession(igvArgs.getSessionFile(), igvArgs.getLocusString());
                            } else {
                                File sf = new File(igvArgs.getSessionFile());
                                if (sf.exists()) {
                                    success = loadSession(sf.getAbsolutePath(), igvArgs.getLocusString());
                                }
                            }
                        } finally {
                            IGV.getInstance().setStatusBarMessage3("");
                        }
                        if (!success) {
                            String genomeId = preferences.getDefaultGenome();
                            //contentPane.getCommandBar().selectGenome(genomeId);
                            try {
                                GenomeManager.getInstance().loadGenomeById(genomeId);
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        }
                    } else if (igvArgs.getDataFileStrings() != null) {

                        // Not an xml file, assume its a list of data files
                        List<String> dataFiles = igvArgs.getDataFileStrings();

                        Collection<String> h = igvArgs.getHttpHeader();
                        if (h != null && !h.isEmpty()) {
                            HttpUtils.getInstance().addHeaders(h, dataFiles);
                        }

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
                        for (int i = 0; i < dataFiles.size(); i++) {

                            String p = dataFiles.get(i).trim();

                            // Decode local file urls??? I don't understand this extra decoding
                            if (URLUtils.isURL(p) && !FileUtils.isRemote(p)) {
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
                                if (URLUtils.isURL(idxP) && !FileUtils.isRemote(idxP)) {
                                    idxP = StringUtils.decodeURL(idxP);       // ???
                                }
                                if (idxP.length() > 0) {
                                    rl.setIndexPath(idxP);
                                }
                            }

                            //Set coverage file, iff one was passed
                            if (coverageFiles != null && i < coverageFiles.length) {
                                String covP = coverageFiles[i];
                                if (URLUtils.isURL(covP) && !FileUtils.isRemote(covP)) {
                                    covP = StringUtils.decodeURL(covP);       // ???
                                }
                                if (covP.length() > 0) {
                                    rl.setCoverage(covP);
                                }
                            }

                            locators.add(rl);
                        }
                        loadTracks(locators);
                    } else if (loadAutosave) {
                        boolean success = false;
                        try {
                            // Get the last autosave and attempt to load
                            File sessionAutosave = SessionAutosaveManager.getMostRecentAutosaveFile().get();
                            success = loadSession(sessionAutosave.getAbsolutePath(), null);
                        } catch (Exception e) {
                            log.error("Failure trying to load most recent autosave file", e);
                        }

                        // Load the default genome if unsuccessful
                        if (!success) {
                            String genomeId = preferences.getDefaultGenome();
                            //contentPane.getCommandBar().selectGenome(genomeId);
                            try {
                                GenomeManager.getInstance().loadGenomeById(genomeId);
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        }
                    }
                }

                UIUtilities.invokeAndWaitOnEventThread(() -> {

                    if (igvArgs.getLocusString() != null) {
                        goToLocus(igvArgs.getLocusString());
                    }

                });

                session.recordHistory();
            }

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


    }

    public static void copySequenceToClipboard(Genome genome, String chr, int start, int end, Strand strand) {
        try {
            IGV.getInstance().getMainFrame().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
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
            getInstance().getMainFrame().setCursor(Cursor.getDefaultCursor());
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

    /**
     * Post an event to this instance's event bus
     * // TODO -- replace the reference to the global event bus with a local one (member if IGV)
     *
     * @param event
     */
    public void postEvent(IGVEvent event) {
        IGVEventBus.getInstance().post(event);
    }


    public void receiveEvent(IGVEvent event) {
        if (event instanceof ViewChange || event instanceof InsertionSelectionEvent) {
            repaint();
        } else if (event instanceof GenomeChangeEvent) {
            repaint();
        } else {
            log.warn("Unknown event type: " + event.getClass());
        }
    }


    public void resetFrames() {
        UIUtilities.invokeOnEventThread(() -> {
                    getMainPanel().headerPanelContainer.createHeaderPanels();
                    for (TrackPanel tp : getTrackPanels()) {
                        tp.createDataPanels();
                    }
                    contentPane.getCommandBar().setGeneListMode(FrameManager.isGeneListMode());
                    revalidateTrackPanels();
                }
        );
    }

    public void revalidateTrackPanels() {
        UIUtilities.invokeOnEventThread(() -> {
            getMainPanel().revalidateTrackPanels();
            repaint(rootPane);
        });
    }


    public void repaintNamePanels() {
        for (TrackPanel tp : getTrackPanels()) {
            tp.getScrollPane().getNamePanel().repaint();
        }
    }


    /**
     * Adjust the height of tracks so that all tracks fit in the available
     * height of the panel. This is not possible in all cases as the
     * minimum height for tracks is respected.
     */
    public void fitTracksToPanel() {
        for (TrackPanel tp : getTrackPanels()) {
            tp.fitTracksToPanel();
        }
        repaint();
    }

    public void repaint() {
        repaint(contentPane);
    }

    public void repaint(Track track) {
        this.repaint(contentPane, List.of(track));
    }

    public void repaint(Collection<? extends Track> tracks) {
        this.repaint(contentPane, tracks);
    }

    private void repaint(final JComponent component) {
        Collection<Track> trackList = new ArrayList<>();
        for (TrackPanel tp : getTrackPanels()) {
            trackList.addAll(visibleTracks(tp.getDataPanelContainer()));
        }
        repaint(component, trackList);
    }

    private boolean isLoading = false;
    private Collection<? extends Track> pending = null;

    private void repaint(final JComponent component, Collection<? extends Track> trackList) {

        if (Globals.isBatch()) {
            // In batch mode everything is done synchronously on the event thread
            UIUtilities.invokeAndWaitOnEventThread(() -> {

                CursorToken token = WaitCursorManager.showWaitCursor();
                try {
                    for (ReferenceFrame frame : FrameManager.getFrames()) {
                        for (Track track : trackList) {
                            if (!track.isReadyToPaint(frame)) {
                                track.load(frame);
                            }
                        }
                    }
                    Autoscaler.autoscale(getAllTracks());
                    checkPanelLayouts();
                    component.paintImmediately(component.getBounds());
                } finally {
                    WaitCursorManager.removeWaitCursor(token);
                    synchronized (IGV.getInstance()) {
                        IGV.getInstance().notifyAll();
                    }
                }
            });

        } else {

            if (isLoading) {
                // Track data is being loaded, do a repaint with existing data and mark this request for future execution
                UIUtilities.invokeOnEventThread(() -> contentPane.repaint());
                pending = trackList;
                return;
            }

            List<CompletableFuture<Void>> futures = new ArrayList<>();

            for (ReferenceFrame frame : FrameManager.getFrames()) {
                for (Track track : trackList) {
                    if (!track.isReadyToPaint(frame)) {
                        futures.add(CompletableFuture.runAsync(() -> track.load(frame), threadExecutor));
                    }
                }
            }

            if (futures.size() == 0) {
                UIUtilities.invokeOnEventThread(() -> {
                    Autoscaler.autoscale(getAllTracks());
                    checkPanelLayouts();
                    component.repaint();
                });
            } else {
                // One or more tracks require loading before repaint.   Load all needed tracks, autoscale if needed, then
                // repaint.  The autoscale step is key, since tracks can be grouped for autoscaling it is necessary that
                // all data is loaded before any track is repainted.  Otherwise tracks be loaded and painted independently.

                final CompletableFuture[] futureArray = futures.toArray(new CompletableFuture[futures.size()]);
                WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                isLoading = true;
                CompletableFuture.allOf(futureArray).thenApplyAsync(future -> {
                    WaitCursorManager.removeWaitCursor(token);
                    // Autoscale as required, check layouts (for scrollbar changes), and repaint.
                    Autoscaler.autoscale(getAllTracks());
                    UIUtilities.invokeOnEventThread(() -> {
                        checkPanelLayouts();
                        component.repaint();
                        isLoading = false;
                        if (pending != null) {
                            Collection<? extends Track> tmp = pending;
                            pending = null;
                            repaint(tmp);
                        }
                    });
                    return null;
                }).exceptionally(ex -> {
                    WaitCursorManager.removeWaitCursor(token);
                    log.error("Error loading track data", ex);
                    isLoading = false;
                    pending = null;
                    return null;
                });
            }
        }
    }

    private void checkPanelLayouts() {
        for (TrackPanel tp : getTrackPanels()) {
            if (tp.isHeightChanged()) {
                UIUtilities.invokeOnEventThread(() -> tp.revalidate());
            }
        }
    }


    public List<Track> visibleTracks(DataPanelContainer dataPanelContainer) {
        return dataPanelContainer.getTrackGroups().stream().
                filter(TrackGroup::isVisible).
                flatMap(trackGroup -> trackGroup.getVisibleTracks().stream()).
                collect(Collectors.toList());
    }

    /**
     * Stops the scheduled autosave task
     */
    public void stopTimedAutosave() {
        sessionAutosaveTimer.cancel();
    }

    // Thread pool for loading data
    private static final ExecutorService threadExecutor = Executors.newFixedThreadPool(5);
}
