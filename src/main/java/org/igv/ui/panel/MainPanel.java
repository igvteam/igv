package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.PreferencesManager;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.track.TrackType;
import org.igv.ui.IGV;
import org.igv.ui.util.UIUtilities;
import org.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
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

        // Enable drag-and-drop for files and URLs
        new DropTarget(this, DnDConstants.ACTION_COPY_OR_MOVE, this, true);

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
        });
    }

    public void removeHeader() {
        remove(headerScrollPane);
        revalidate();
    }

    public void restoreHeader() {
        add(headerScrollPane, BorderLayout.NORTH);
        revalidate();
    }


    @Override
    public void doLayout() {
        super.doLayout();
        applicationHeaderPanel.doLayout();
        for (TrackPanel tp : getTrackPanels()) {
            tp.getScrollPane().doLayout();
        }
    }

    @Override
    public void setBackground(Color color) {
        super.setBackground(color);
        if (headerPanelContainer != null) {
            applicationHeaderPanel.setBackground(color);
            nameHeaderPanel.setBackground(color);
            attributeHeaderPanel.setBackground(color);
            headerPanelContainer.setBackground(color);
            nameHeaderPanel.setBackground(color);
            attributeHeaderPanel.setBackground(color);
            for (TrackPanel tsp : getTrackPanels()) {
                tsp.setBackground(color);
            }
        }

    }

    private void initComponents() {

        setPreferredSize(new java.awt.Dimension(1021, 510));
        setLayout(new java.awt.BorderLayout());

        nameHeaderPanel = new NameHeaderPanel();
        nameHeaderPanel.setBackground(new java.awt.Color(255, 255, 255));
        nameHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        nameHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));

        attributeHeaderPanel = new AttributeHeaderPanel();
        attributeHeaderPanel.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        attributeHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        attributeHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));


        headerPanelContainer = new HeaderPanelContainer();
        headerScrollPane = new JScrollPane();
        headerScrollPane.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(102, 102, 102)));
        headerScrollPane.setForeground(new java.awt.Color(153, 153, 153));
        headerScrollPane.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        headerScrollPane.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        headerScrollPane.setPreferredSize(new java.awt.Dimension(1021, 130));
        add(headerScrollPane, java.awt.BorderLayout.NORTH);

        applicationHeaderPanel = new IGVPanel(this);

        // Add spacer panel to align with drag handle in TrackPanel
        JPanel dragHandleSpacer = new JPanel();
        dragHandleSpacer.setBackground(new java.awt.Color(255, 255, 255));
        dragHandleSpacer.setPreferredSize(new java.awt.Dimension(DragHandlePanel.DRAG_HANDLE_WIDTH, 0));
        applicationHeaderPanel.add(dragHandleSpacer);

        applicationHeaderPanel.add(nameHeaderPanel);
        applicationHeaderPanel.add(attributeHeaderPanel);
        applicationHeaderPanel.add(headerPanelContainer);
        headerScrollPane.setViewportView(applicationHeaderPanel);


        // Custom panel that implements Scrollable to prevent viewport from stretching it
        trackPanelContainer = new ScrollableTrackContainer();

        trackPanelScrollPane = new JScrollPane(trackPanelContainer);
        trackPanelScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        trackPanelScrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
        trackPanelScrollPane.setBorder(null);
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

    public synchronized TrackPanelScrollPane addTrackPanel(Track track) {

        final TrackPanel trackPanel = new TrackPanel(track.getName(), this);

        TrackType trackType = track.getTrackType();

        trackPanel.addTrack(track);
        final TrackPanelScrollPane sp = new TrackPanelScrollPane();
        track.setViewport(sp);

        Runnable runnable = () -> {
            sp.setViewportView(trackPanel);
            if (trackType == TrackType.SEQUENCE) {
                trackPanelContainer.add(sp, 0);
            } else {
                int position = 0;
                Component[] components = trackPanelContainer.getComponents();
                if (components.length > 0) {
                    Component c = components[0];
                    if (c instanceof TrackPanelScrollPane) {
                        TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;
                        if (!tsp.getTrackPanel().getTracks().isEmpty()) {
                            Track firstTrack = tsp.getTrackPanel().getTracks().get(0);
                            if (firstTrack.getTrackType() == TrackType.SEQUENCE) {
                                position = 1;
                            }
                        }
                    }
                }
                trackPanelContainer.add(sp, position);
            }
        };

        UIUtilities.invokeAndWaitOnEventThread(runnable);

        return sp;
    }

    public void clearTrackPanels() {
        trackPanelContainer.removeAll();
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

    public void reorderPanels(java.util.List<String> names) {

        Map<String, TrackPanelScrollPane> panes = new HashMap();
        for (Component c : trackPanelContainer.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;
                panes.put(tsp.getTrackPanelName(), tsp);
            }
        }

        trackPanelContainer.removeAll();
        for (String name : names) {
            trackPanelContainer.add(panes.get(name));
        }
        trackPanelContainer.revalidate();
        trackPanelContainer.repaint();
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
                TrackNamePanel.removeDropListenerFor(panel.getNamePanel());
            }

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
            TrackNamePanel.removeDropListenerFor(sp.getNamePanel());
            trackPanelContainer.revalidate();
        }
    }


    public boolean panelIsRemovable(TrackPanel trackPanel) {
        return true;
    }

    public void updatePanelDimensions() {
        Insets insets = applicationHeaderPanel.getInsets();
        namePanelX = insets.left;
        attributePanelX = namePanelX + namePanelWidth + hgap;
        attributePanelWidth = calculateAttributeWidth();
        dataPanelX = attributePanelX + attributePanelWidth + hgap;
        dataPanelWidth = applicationHeaderPanel.getWidth() - insets.right - dataPanelX;
    }

    public int calculateAttributeWidth() {

        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY)) {
            return 0;
        }
        Collection<String> attributeKeys = AttributeManager.getInstance().getVisibleAttributes();
        int attributeCount = attributeKeys.size();
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

    public ScrollableTrackContainer getTrackPanelContainer() {
        return trackPanelContainer;
    }

    @Override
    protected void paintComponent(Graphics g) {
        if (Globals.isDarkMode()) {
            setBackground(UIManager.getColor("Panel.background"));
        }
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
        Rectangle r = trackPanelContainer.getBounds();
        g.translate(0, r.y);

        // Get the components of the center pane and sort by Y position.
        Component[] components = trackPanelContainer.getComponents();
        Arrays.sort(components, Comparator.comparingInt(Component::getY));

        int dy = components[0].getY();
        for (Component c : components) {

            Graphics2D g2d = (Graphics2D) g.create();
            g2d.translate(0, dy);

            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                //Skip if panel has no tracks
                if (tsp.getTrackPanel().getTracks().size() == 0) {
                    continue;
                }

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

//        //super.paintBorder(g);

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

            // Load dropped files and URLs as tracks
            List<ResourceLocator> locators = new ArrayList<>();

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
}
