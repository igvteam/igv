package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.PreferencesManager;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.track.DataType;
import org.igv.ui.IGV;
import org.igv.ui.util.UIUtilities;
import org.igv.util.LongRunningTask;
import org.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.ContainerEvent;
import java.awt.event.ContainerListener;
import java.io.File;
import java.net.URI;
import java.util.*;
import java.util.List;

import static org.igv.prefs.Constants.*;

/**
 * @author jrobinso
 * @date Sep 10, 2010
 */
public class MainPanel extends JPanel implements Paintable, DropTargetListener {

    private static Logger log = LogManager.getLogger(MainPanel.class);

    IGV igv;

    // private static final int DEFAULT_NAME_PANEL_WIDTH = 160;

    private int namePanelX;
    private int namePanelWidth = PreferencesManager.getPreferences().getAsInt(NAME_PANEL_WIDTH);
    private int attributePanelX;
    private int attributePanelWidth;
    private int dataPanelX;
    private int dataPanelWidth;

    public IGVPanel applicationHeaderPanel;
    public HeaderPanelContainer headerPanelContainer;
    private ScrollableTrackContainer trackPanelContainer;
    private JScrollPane trackPanelScrollPane;
    private NameHeaderPanel nameHeaderPanel;
    private AttributeHeaderPanel attributeHeaderPanel;

    private int hgap = 5;
    private JScrollPane headerScrollPane;


    public MainPanel(IGV igv) {
        this.igv = igv;

        initComponents();

        // Enable drag-and-drop for files and URLs.
        // DropTarget events don't bubble from child to parent in AWT/Swing, so we
        // install MainPanel's own DropTarget first, then recursively install
        // forwarding DropTargets on all descendants that don't already have one.
        new DropTarget(this, DnDConstants.ACTION_COPY_OR_MOVE, this, true);
        installDropTargetRecursively(this);

        //Load IGV logo
//        try {
//            BufferedImage logo = ImageIO.read(getClass().getResource("resources/IGV_64.png"));
//            JLabel picLabel = new JLabel(new ImageIcon(logo));
//            picLabel.setVerticalAlignment(SwingConstants.CENTER);
//            nameHeaderPanel.add(picLabel);
//        } catch (IOException e) {
//            //pass
//        }

        addComponentListener(new ComponentListener() {

            public void componentResized(ComponentEvent componentEvent) {
                revalidateTrackPanels();
                igv.repaint();
            }

            public void componentMoved(ComponentEvent componentEvent) {
                //To change body of implemented methods use File | Settings | File Templates.
            }

            public void componentShown(ComponentEvent componentEvent) {
                //To change body of implemented methods use File | Settings | File Templates.
            }

            public void componentHidden(ComponentEvent componentEvent) {
                //To change body of implemented methods use File | Settings | File Templates.
            }
        });
    }


    public void collapseNamePanel() {
        namePanelWidth = 0;
        revalidateTrackPanels();
    }

    public void expandNamePanel() {
        namePanelWidth = PreferencesManager.getPreferences().getAsInt(NAME_PANEL_WIDTH);
        revalidateTrackPanels();
    }

    public void setNamePanelWidth(int width) {
        this.namePanelWidth = width;
        revalidateTrackPanels();
    }

    public void revalidateTrackPanels() {
        updatePanelDimensions();
        UIUtilities.invokeOnEventThread(() -> {
            this.applicationHeaderPanel.invalidate();
            for (TrackPanel tp : this.getTrackPanels()) {
                tp.invalidate();
            }
            this.invalidate(); // this should not be neccessary, but is harmless
            this.validate();
            this.repaint();  // Repaint to update divider lines (horizontal separator + track area vertical dividers)
        });
    }

    public void removeHeader() {
        remove(headerScrollPane);
        trackPanelScrollPane.setBorder(null);
        revalidate();
    }

    public void restoreHeader() {
        add(headerScrollPane, BorderLayout.NORTH);
        trackPanelScrollPane.setBorder(createHeaderSeparatorBorder());
        revalidate();
    }

    /**
     * Creates a border with a 1px top line to visually separate the header from the track area.
     */
    private static javax.swing.border.Border createHeaderSeparatorBorder() {
        Color dividerColor = Globals.isDarkMode() ? Color.GRAY : Color.LIGHT_GRAY;
        return BorderFactory.createMatteBorder(1, 0, 0, 0, dividerColor);
    }


    @Override
    public void doLayout() {
        super.doLayout();
        applicationHeaderPanel.doLayout();
        for (TrackPanel tp : getTrackPanels()) {
            tp.getScrollPane().doLayout();
        }
    }


    private void initComponents() {

        setPreferredSize(new java.awt.Dimension(1021, 510));
        setLayout(new java.awt.BorderLayout());

        nameHeaderPanel = new NameHeaderPanel();
        nameHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        nameHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));

        attributeHeaderPanel = new AttributeHeaderPanel();
        attributeHeaderPanel.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        attributeHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        attributeHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));


        headerPanelContainer = new HeaderPanelContainer();
        headerScrollPane = new JScrollPane();
        headerScrollPane.setBorder(null);
        // headerScrollPane.setForeground(new java.awt.Color(153, 153, 153));
        headerScrollPane.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        headerScrollPane.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        headerScrollPane.setPreferredSize(new java.awt.Dimension(1021, 130));
        add(headerScrollPane, java.awt.BorderLayout.NORTH);

        applicationHeaderPanel = new IGVPanel(this);

        // Add spacer panel to align with drag handle in TrackPanel
        JPanel dragHandleSpacer = new JPanel();
        //dragHandleSpacer.setBackground(new java.awt.Color(255, 255, 255));
        dragHandleSpacer.setPreferredSize(new java.awt.Dimension(DragHandlePanel.DRAG_HANDLE_WIDTH, 0));
        applicationHeaderPanel.add(dragHandleSpacer);

        applicationHeaderPanel.add(nameHeaderPanel);
        applicationHeaderPanel.add(attributeHeaderPanel);
        applicationHeaderPanel.add(headerPanelContainer);
        headerScrollPane.setViewportView(applicationHeaderPanel);


        // Custom panel that implements Scrollable to prevent viewport from stretching it
        trackPanelContainer = new ScrollableTrackContainer(this);

        trackPanelScrollPane = new JScrollPane(trackPanelContainer);
        trackPanelScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        trackPanelScrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
        trackPanelScrollPane.setBorder(createHeaderSeparatorBorder());
        trackPanelScrollPane.getVerticalScrollBar().setUnitIncrement(16);
        trackPanelScrollPane.getVerticalScrollBar().addAdjustmentListener(e -> {
            if (!e.getValueIsAdjusting()) {
                Rectangle viewRect = trackPanelScrollPane.getViewport().getViewRect();
                for (Component c : trackPanelContainer.getComponents()) {
                    if (c instanceof TrackPanelScrollPane) {
                        TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;
                        if (c.getBounds().intersects(viewRect)) {
                            tsp.trackPanel.getNamePanel().repaint();
                        }
                    }
                }
            }
        });
        add(trackPanelScrollPane, BorderLayout.CENTER);

        setBackground(PreferencesManager.getPreferences().getAsColor(BACKGROUND_COLOR));


    }

    /**
     * Add a track panel at the position determined by the track's order property.
     * Tracks are inserted to maintain ascending order by the order property.
     *
     * @param track the track to add
     * @return the TrackPanelScrollPane containing the track
     */
    public synchronized TrackPanelScrollPane addTrackPanel(Track track) {

        final TrackPanel trackPanel = new TrackPanel(track.getName(), this);

        trackPanel.addTrack(track);
        final TrackPanelScrollPane sp = new TrackPanelScrollPane();
        track.setViewport(sp);

        Runnable runnable = () -> {
            sp.setViewportView(trackPanel);
            long trackOrder = track.getOrder();
            if (trackOrder == 0) {
                track.setOrder(getTrackPanels().size());
            }
            int insertPosition = findInsertPosition(track.getOrder());
            trackPanelContainer.add(sp, insertPosition);
        };

        UIUtilities.invokeAndWaitOnEventThread(runnable);

        rebuildDividers();

        return sp;
    }

    /**
     * Find the correct insertion position for a track with the given order value.
     * Returns the index where the track should be inserted to maintain ascending order.
     *
     * @param order the order value of the track to insert
     * @return the index at which to insert the track panel
     */
    private int findInsertPosition(long order) {
        Component[] components = trackPanelContainer.getComponents();
        for (int i = 0; i < components.length; i++) {
            Component c = components[i];
            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;
                List<Track> tracks = tsp.getTrackPanel().getTracks();
                if (!tracks.isEmpty()) {
                    Track existingTrack = tracks.get(0);
                    if (existingTrack.getOrder() > order) {
                        return i;
                    }
                }
            }
        }
        // If no track has a higher order, insert at the end
        return components.length;
    }

    /**
     * Updates the order property of a track that was moved via drag & drop.
     * Calculates a new order value that places it between its new neighbors,
     * preserving any pinned tracks with extreme order values.
     *
     * @param movedPanel the track panel that was moved
     */
    public void updateMovedTrackOrder(TrackPanel movedPanel) {
        List<TrackPanel> trackPanels = getTrackPanels();
        int newIndex = trackPanels.indexOf(movedPanel);
        if (newIndex < 0) return;

        Track movedTrack = movedPanel.getTrack();
        if (movedTrack == null) return;

        long prevOrder = Globals.JS_MIN_SAFE_INTEGER;
        long nextOrder = Globals.JS_MAX_SAFE_INTEGER;

        // Get the order of the previous track (if any)
        if (newIndex > 0) {
            Track prevTrack = trackPanels.get(newIndex - 1).getTrack();
            if (prevTrack != null) {
                prevOrder = prevTrack.getOrder();
            }
        }

        // Get the order of the next track (if any)
        if (newIndex < trackPanels.size() - 1) {
            Track nextTrack = trackPanels.get(newIndex + 1).getTrack();
            if (nextTrack != null) {
                nextOrder = nextTrack.getOrder();
            }
        }

        // Calculate new order as midpoint between neighbors
        long newOrder;
        if (prevOrder == Globals.JS_MIN_SAFE_INTEGER && nextOrder == Globals.JS_MAX_SAFE_INTEGER) {
            // Only track, use 0
            newOrder = 0;
        } else if (prevOrder == Globals.JS_MIN_SAFE_INTEGER) {
            // No previous track, place before next
            newOrder = nextOrder - 1;
        } else if (nextOrder == Globals.JS_MAX_SAFE_INTEGER) {
            // No next track, place after previous
            newOrder = prevOrder + 1;
        } else {
            // Between two tracks, use midpoint
            newOrder = prevOrder + (nextOrder - prevOrder) / 2;
            // If midpoint equals prevOrder (no room), just use prevOrder + 1
            if (newOrder == prevOrder) {
                newOrder = prevOrder + 1;
            }
        }

        movedTrack.setOrder(newOrder);
    }

    /**
     * Get the index of the track panel containing the specified track.
     *
     * @param track the track to find
     * @return the index of the track panel, or -1 if not found
     */
    public int getTrackPanelIndex(Track track) {
        Component[] components = trackPanelContainer.getComponents();
        for (int i = 0; i < components.length; i++) {
            Component c = components[i];
            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;
                if (tsp.getTrackPanel().containsTrack(track)) {
                    return i;
                }
            }
        }
        return -1;
    }

    public void clearTrackPanels() {
        trackPanelContainer.removeAll();
        rebuildDividers();
    }

    /**
     * Return an ordered list of TrackPanels.  This method is provided primarily for storing sessions, where
     * TrackPanels need to be stored in proper order
     *
     * @return
     */
    public java.util.List<TrackPanel> getTrackPanels() {
        ArrayList<TrackPanel> panels = new ArrayList<TrackPanel>();
        for (Component c : trackPanelContainer.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                panels.add(((TrackPanelScrollPane) c).getTrackPanel());
            }
        }
        return panels;
    }

    /**
     * Reorder track panels using direct scroll pane references.
     * This avoids any name-matching issues that could cause panels to be lost.
     *
     * @param orderedPanes the scroll panes in the desired order
     */
    public void reorderPanels(java.util.List<TrackPanelScrollPane> orderedPanes) {
        trackPanelContainer.removeAll();
        for (TrackPanelScrollPane pane : orderedPanes) {
            if (pane != null) {
                trackPanelContainer.add(pane);
            }
        }
        rebuildDividers();
    }

    /**
     * Reorder track panels by name. Used by {@link ReorderPanelsDialog}.
     * Delegates to {@link #reorderPanels(List)} with resolved scroll pane references.
     *
     * @param names the panel names in the desired order
     */
    public void reorderPanelsByName(java.util.List<String> names) {
        Map<String, TrackPanelScrollPane> panesByName = new HashMap<>();
        for (Component c : trackPanelContainer.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;
                panesByName.put(tsp.getTrackPanelName(), tsp);
            }
        }

        java.util.List<TrackPanelScrollPane> orderedPanes = new ArrayList<>();
        for (String name : names) {
            TrackPanelScrollPane pane = panesByName.get(name);
            if (pane != null) {
                orderedPanes.add(pane);
            }
        }
        reorderPanels(orderedPanes);
    }

    public void removeEmptyDataPanels() {
        List<TrackPanelScrollPane> emptyPanels = new ArrayList();
        for (TrackPanel tp : getTrackPanels()) {
            if (tp.getTracks().isEmpty()) {
                emptyPanels.add(tp.getScrollPane());
            }
        }
        for (TrackPanelScrollPane panel : emptyPanels) {
            if (panel != null) {
                trackPanelContainer.remove(panel);
            }
        }
        if (!emptyPanels.isEmpty()) {
            rebuildDividers();
        }
    }

    public void removeDataPanel(String name) {

        for (TrackPanel tp : getTrackPanels()) {
            if (name.equals(tp.getName())) {
                removeTrackPanel(tp);
                return;
            }
        }
    }

    public void removeTrackPanel(TrackPanel trackPanel) {
        TrackPanelScrollPane sp = trackPanel.getScrollPane();
        if (sp != null) {
            trackPanelContainer.remove(sp);
            rebuildDividers();
            trackPanelContainer.revalidate();
        }
    }


    public boolean panelIsRemovable(TrackPanel trackPanel) {
        return true;
    }

    /**
     * Rebuild dividers between TrackPanelScrollPanes. Removes all existing
     * {@link TrackPanelDivider} instances from the container and inserts a new
     * divider after each TrackPanelScrollPane (including the last one, so that
     * the last track can also be resized).
     */
    private void rebuildDividers() {
        // Collect the current TrackPanelScrollPanes in order
        List<TrackPanelScrollPane> panes = new ArrayList<>();
        for (Component c : trackPanelContainer.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                panes.add((TrackPanelScrollPane) c);
            }
        }

        // Re-add panes, each followed by a divider
        trackPanelContainer.removeAll();
        for (int i = 0; i < panes.size(); i++) {
            TrackPanelScrollPane above = panes.get(i);
            TrackPanelScrollPane below = (i + 1 < panes.size()) ? panes.get(i + 1) : null;
            trackPanelContainer.add(above);
            trackPanelContainer.add(new TrackPanelDivider(above, below));
        }

        trackPanelContainer.revalidate();
        trackPanelContainer.repaint();
    }

    public void updatePanelDimensions() {
        Insets insets = applicationHeaderPanel.getInsets();
        namePanelX = insets.left;
        attributePanelWidth = calculateAttributeWidth();
        if (attributePanelWidth > 0) {
            attributePanelX = namePanelX + namePanelWidth + hgap;
            dataPanelX = attributePanelX + attributePanelWidth + hgap;
        } else {
            attributePanelX = namePanelX + namePanelWidth + hgap;
            dataPanelX = attributePanelX;
        }
        dataPanelWidth = applicationHeaderPanel.getWidth() - insets.right - dataPanelX;
    }

    public int calculateAttributeWidth() {

        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY)) {
            return 0;
        }
        Collection<String> attributeKeys = AttributeManager.getInstance().getVisibleAttributes();
        int attributeCount = attributeKeys.size();
        if (attributeCount == 0) {
            return 0;
        }
        int packWidth = (attributeCount) * (AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH +
                AttributeHeaderPanel.COLUMN_BORDER_WIDTH) + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
        return packWidth;
    }

    public boolean isExpanded() {
        return namePanelWidth > 0;
    }

    public int getAttributePanelWidth() {
        return attributePanelWidth;
    }

    public int getNamePanelX() {
        return namePanelX;
    }

    public int getNamePanelWidth() {
        return namePanelWidth;
    }

    public int getAttributePanelX() {
        return attributePanelX;
    }


    public int getDataPanelX() {
        return dataPanelX;
    }


    public int getDataPanelWidth() {
        return dataPanelWidth;
    }

    /**
     * Returns the total left offset for panels, which includes the drag handle width
     * plus the selection panel width if selection panels are visible.
     */
    public int getLeftOffset() {
        int offset = DragHandlePanel.DRAG_HANDLE_WIDTH;
        if (TrackSelectionPanel.isSelectionModeActive()) {
            offset += TrackSelectionPanel.SELECTION_PANEL_WIDTH;
        }
        return offset;
    }

    public ScrollableTrackContainer getTrackPanelContainer() {
        return trackPanelContainer;
    }



    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        // Header
        int width = applicationHeaderPanel.getWidth();
        int height = applicationHeaderPanel.getHeight();

        Graphics2D headerGraphics = (Graphics2D) g.create();
        Rectangle headerRect = new Rectangle(0, 0, width, height);
        applicationHeaderPanel.paintOffscreen(headerGraphics, headerRect, batch);
        headerGraphics.dispose();

        // Now loop through track panels
        g.translate(0, height);

        // Get the components of the center pane and sort by Y position.
        Component[] components = trackPanelContainer.getComponents();
        Arrays.sort(components, Comparator.comparingInt(Component::getY));

        int dy = components[0].getY();
        for (Component c : components) {

            Graphics2D g2d = (Graphics2D) g.create();
            g2d.translate(0, dy);

            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                int panelHeight = tsp.getSnapshotHeight(batch);

                Rectangle tspRect = new Rectangle(0, 0, tsp.getWidth(), panelHeight);

                g2d.setClip(tspRect);
                tsp.paintOffscreen(g2d, tspRect, batch);
                dy += tspRect.height;

            } else {
                g2d.setClip(new Rectangle(0, 0, c.getWidth(), c.getHeight()));
                c.paint(g2d);
                dy += c.getHeight();
            }

            g2d.dispose();

        }
    }

    /**
     * Return the image height required to paint this component with current options.  This is used to size bitmap
     * images for offscreen drawing.
     *
     * @return
     */
    @Override
    public int getSnapshotHeight(boolean batch) {

        if (batch) {
            int height = applicationHeaderPanel.getHeight();

            for (Component c : trackPanelContainer.getComponents()) {

                if (c instanceof TrackPanelScrollPane) {

                    TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                    //Skip if panel has no tracks
                    if (tsp.getTrackPanel().getTracks().size() == 0) {
                        continue;
                    }

                    height += tsp.getSnapshotHeight(batch);

                } else {
                    height += c.getHeight();
                }

            }
            return height;
        } else {
            return getHeight();
        }
    }

    public void repaintHeaderPanels() {
        headerPanelContainer.repaint();
    }

    // DropTargetListener implementation

    @Override
    public void dragEnter(DropTargetDragEvent dtde) {
        if (isDragAcceptable(dtde)) {
            dtde.acceptDrag(DnDConstants.ACTION_COPY);
        } else {
            dtde.rejectDrag();
        }
    }

    @Override
    public void dragOver(DropTargetDragEvent dtde) {
        // No action needed
    }

    @Override
    public void dropActionChanged(DropTargetDragEvent dtde) {
        // No action needed
    }

    @Override
    public void dragExit(DropTargetEvent dte) {
        // No action needed
    }

    @Override
    public void drop(DropTargetDropEvent dtde) {
        try {
            dtde.acceptDrop(DnDConstants.ACTION_COPY);
            Transferable transferable = dtde.getTransferable();

            List<File> droppedFiles = new ArrayList<>();
            List<String> droppedUrls = new ArrayList<>();

            // Try to get files
            if (transferable.isDataFlavorSupported(DataFlavor.javaFileListFlavor)) {
                @SuppressWarnings("unchecked")
                List<File> files = (List<File>) transferable.getTransferData(DataFlavor.javaFileListFlavor);
                droppedFiles.addAll(files);
            }

            // Try to get URLs/URIs (for URL drops, including from browsers)
            if (transferable.isDataFlavorSupported(DataFlavor.stringFlavor)) {
                String data = (String) transferable.getTransferData(DataFlavor.stringFlavor);
                if (data != null && !data.trim().isEmpty()) {
                    // Parse potential URLs (one per line)
                    String[] lines = data.split("\\r?\\n");
                    for (String line : lines) {
                        line = line.trim();
                        if (isUrl(line)) {
                            droppedUrls.add(line);
                        }
                    }
                }
            }

            // Try URI list flavor (common on Linux)
            DataFlavor uriListFlavor = new DataFlavor("text/uri-list;class=java.lang.String");
            if (transferable.isDataFlavorSupported(uriListFlavor)) {
                String data = (String) transferable.getTransferData(uriListFlavor);
                if (data != null) {
                    String[] uris = data.split("\\r?\\n");
                    for (String uriStr : uris) {
                        uriStr = uriStr.trim();
                        if (!uriStr.isEmpty() && !uriStr.startsWith("#")) {
                            try {
                                URI uri = new URI(uriStr);
                                if ("file".equals(uri.getScheme())) {
                                    droppedFiles.add(new File(uri));
                                } else {
                                    droppedUrls.add(uriStr);
                                }
                            } catch (Exception e) {
                                // Ignore invalid URIs
                            }
                        }
                    }
                }
            }

            // Load dropped files and URLs as tracks or sessions
            List<ResourceLocator> locators = new ArrayList<>();

            // First check for a session
            String sessionPath = null;
            if (droppedFiles.size() == 1 && droppedFiles.get(0).getName().endsWith(".xml")) {
                // If it's a single XML file, treat it as a session file
                sessionPath = droppedFiles.get(0).getAbsolutePath();
            } else if (droppedUrls.size() == 1 && droppedUrls.get(0).endsWith(".xml")) {
                // If it's a single URL ending in .xml, treat it as a session URL
                sessionPath = droppedUrls.get(0);
            }
            if (sessionPath != null) {
                final String sp = sessionPath;
                LongRunningTask.submit(() -> this.igv.loadSession(sp, null));
                return;
            }

            // Add file locators
            if (!droppedFiles.isEmpty()) {
                locators.addAll(ResourceLocator.getLocators(droppedFiles));
            }

            // Add URL locators
            for (String url : droppedUrls) {
                locators.add(new ResourceLocator(url));
            }

            if (!locators.isEmpty()) {
                igv.loadTracks(locators);
            }

            dtde.dropComplete(true);
        } catch (Exception e) {
            dtde.dropComplete(false);
        }
    }

    private boolean isDragAcceptable(DropTargetDragEvent dtde) {
        return dtde.isDataFlavorSupported(DataFlavor.javaFileListFlavor) ||
                dtde.isDataFlavorSupported(DataFlavor.stringFlavor);
    }

    private boolean isUrl(String str) {
        if (str == null) return false;
        String lower = str.toLowerCase();
        return lower.startsWith("http://") ||
                lower.startsWith("https://") ||
                lower.startsWith("ftp://") ||
                lower.startsWith("s3://") ||
                lower.startsWith("gs://");
    }

    /**
     * Recursively install a DropTarget on the given component and all its
     * descendants so that drop events anywhere in the UI are handled.
     * <p>
     * Components that already have their own DropTarget (e.g. TrackPanel for reorder D&D,
     * HeaderPanel for frame reorder D&D) are left untouched — their listeners already
     * delegate unhandled flavors (like file drops) to {@link #drop}.
     * <p>
     * For components without a DropTarget, a forwarding listener is installed that
     * locates the nearest ancestor with a DropTarget and dispatches the event there.
     * This preserves the existing D&D behaviour (e.g. track reorder) while also making
     * file drops work everywhere.
     * <p>
     * A ContainerListener is added to every Container so that dynamically added
     * children (e.g. new track panels) are handled automatically.
     */
    private void installDropTargetRecursively(Component comp) {
        if (comp.getDropTarget() == null) {
            // Install a forwarding DropTarget that bubbles the event up to the nearest
            // ancestor that has its own DropTarget (ultimately reaching MainPanel).
            DropTargetListener forwarder = new DropTargetAdapter() {
                @Override
                public void drop(DropTargetDropEvent dtde) {
                    DropTarget ancestorDt = findAncestorDropTarget(comp);
                    if (ancestorDt != null) {
                        ancestorDt.drop(dtde);
                    }
                }
            };
            new DropTarget(comp, DnDConstants.ACTION_COPY_OR_MOVE, forwarder, true);
        }
        if (comp instanceof Container) {
            Container container = (Container) comp;
            for (Component child : container.getComponents()) {
                installDropTargetRecursively(child);
            }
            container.addContainerListener(new ContainerListener() {
                @Override
                public void componentAdded(ContainerEvent e) {
                    installDropTargetRecursively(e.getChild());
                }

                @Override
                public void componentRemoved(ContainerEvent e) {
                    // No action needed on removal
                }
            });
        }
    }

    /**
     * Walk up the component hierarchy and return the DropTarget of the first
     * ancestor that has one.  Returns {@code null} if none is found.
     */
    private static DropTarget findAncestorDropTarget(Component comp) {
        Container parent = comp.getParent();
        while (parent != null) {
            DropTarget dt = parent.getDropTarget();
            if (dt != null) {
                return dt;
            }
            parent = parent.getParent();
        }
        return null;
    }
}
