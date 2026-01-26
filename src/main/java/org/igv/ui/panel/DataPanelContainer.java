package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.feature.genome.GenomeManager;
import org.igv.feature.genome.load.HubGenomeLoader;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.session.SessionReader;
import org.igv.track.Track;
import org.igv.track.TrackGroup;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;
import org.igv.ui.MessageCollection;
import org.igv.ui.util.MessageUtils;
import org.igv.util.LongRunningTask;
import org.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 8, 2010
 */
public class DataPanelContainer extends TrackPanelComponent implements Paintable {

    static final int default_hgap = 6;
    private static Logger log = LogManager.getLogger(DataPanelContainer.class);

    TrackPanel parent;

    public DataPanelContainer(TrackPanel trackPanel) {
        super(trackPanel);
        DropTarget target = new DropTarget(this, new FileDropTargetListener(trackPanel));
        setDropTarget(target);
        target.setActive(true);
        this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
        this.parent = trackPanel;
        createDataPanels();
    }

    public void createDataPanels() {
        removeAll();
        int hgap = default_hgap;
        if (FrameManager.getFrames().size() > 10) {
            hgap = 1 + 20 / FrameManager.getFrames().size();
        }
        boolean first = true;
        for (ReferenceFrame f : FrameManager.getFrames()) {
            if (f.isVisible()) {
                if (!first) {
                    this.add(Box.createRigidArea(new Dimension(hgap, 0)));
                }
                DataPanel dp = new DataPanel(f, this);
                add(dp);
                first = false;
            }
        }
        invalidate();
    }

    @Override
    public void setBackground(Color color) {
        super.setBackground(color);
        for (Component c : this.getComponents()) {
            if (c instanceof DataPanel) {
                c.setBackground(color);
            }
        }
    }


    public Collection<TrackGroup> getTrackGroups() {
        TrackPanel dataTrackView = (TrackPanel) getParent();
        return dataTrackView.getGroups();
    }

    public int getVisibleHeight() {
        TrackPanel dataTrackView = (TrackPanel) getParent();
        return dataTrackView.getVisibleRect().height;
    }

    public void setCurrentTool(RegionOfInterestTool regionOfInterestTool) {
        for (Component c : this.getComponents()) {
            if (c instanceof DataPanel) {
                ((DataPanel) c).setCurrentTool(regionOfInterestTool);
            }
        }
    }

    /**
     * Paint to an offscreen graphic, e.g. a graphic for an image or svg file.
     *
     * @param g -- graphics context, translated as neccessary to datapanel origin
     * @param rect  -- Rectangle in which to draw datapanel container
     */
    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        // Get the components of the sort by X position.
        Component[] components = getComponents();
        Arrays.sort(components, Comparator.comparingInt(Component::getX));

        for (Component c : this.getComponents()) {
            if (c instanceof DataPanel) {
                Graphics2D g2d = (Graphics2D) g.create();
                g2d.translate(c.getX(), 0);
                Rectangle panelRect = new Rectangle(0, rect.y, c.getWidth(), rect.height);
                ((DataPanel) c).paintOffscreen(g2d, panelRect, batch);
            }
        }
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
    }

    @Override
    protected void paintChildren(Graphics g) {

        super.paintChildren(g);
        if (IGV.getInstance().isRulerEnabled()) {
            int start = MouseInfo.getPointerInfo().getLocation().x - getLocationOnScreen().x;
            g.setColor(Color.BLACK);
            g.drawLine(start, 0, start, getHeight());

            ReferenceFrame frame = FrameManager.getDefaultFrame();
            boolean allChrMode = frame.getChrName().equals(Globals.CHR_ALL);

            if (!FrameManager.isGeneListMode() && !allChrMode) {
                int y = MouseInfo.getPointerInfo().getLocation().y - getLocationOnScreen().y;
                int pos = (int) frame.getChromosomePosition(start) + 1;

                g.setFont(FontManager.getDefaultFont());
                g.drawString(Globals.DECIMAL_FORMAT.format((double) pos), start + 10, y + 30);
            }
        }
    }


    private class FileDropTargetListener implements DropTargetListener {

        private TrackPanel panel;

        public FileDropTargetListener(TrackPanel dataPanel) {
            panel = dataPanel;
        }

        public void dragEnter(DropTargetDragEvent event) {

            if (!isDragAcceptable(event)) {
                event.rejectDrag();
                return;
            }
        }

        public void dragExit(DropTargetEvent event) {
        }

        public void dragOver(DropTargetDragEvent event) {
            // you can provide visual feedback here
        }

        public void dropActionChanged(DropTargetDragEvent event) {
            if (!isDragAcceptable(event)) {
                event.rejectDrag();
                return;
            }
        }

        public void drop(DropTargetDropEvent event) {
            if (!isDropAcceptable(event)) {
                event.rejectDrop();
                return;
            }

            event.acceptDrop(DnDConstants.ACTION_COPY);

            Transferable transferable = event.getTransferable();
            DataFlavor[] flavors = transferable.getTransferDataFlavors();

            MessageCollection messages = new MessageCollection();
            final IGV igv = IGV.getInstance();
            try {
                if (transferable.isDataFlavorSupported(DataFlavor.javaFileListFlavor)) {
                    List<File> files = (List<File>) transferable.getTransferData(DataFlavor.javaFileListFlavor);
                    if (files != null && files.size() > 0) {
                        // Check if this is a single session file
                        if (files.size() == 1 && SessionReader.isSessionFile(files.get(0).getAbsolutePath())) {
                            String sessionPath = files.get(0).getAbsolutePath();
                            LongRunningTask.submit(() -> igv.loadSession(sessionPath, null));
                        } else {
                            List<ResourceLocator> locators = ResourceLocator.getLocators(files);
                            igv.loadTracks(locators, panel);
                        }
                    }
                } else if (transferable.isDataFlavorSupported(DataFlavor.stringFlavor)) {
                    String obj = transferable.getTransferData(DataFlavor.stringFlavor).toString();

                    // Check for genomes, sessions, etc first
                    if (HubGenomeLoader.isHubURL(obj)) {
                        LongRunningTask.submit(() -> {
                            try {
                                GenomeManager.getInstance().loadGenome(obj);
                            } catch (IOException e) {
                                MessageUtils.showMessage("Error loading track hub: " + e.getMessage());
                            }
                        });
                    } else if (SessionReader.isSessionFile(obj)) {
                        LongRunningTask.submit(() -> igv.loadSession(obj, null));
                    } else {
                        igv.loadTracks(Collections.singletonList(new ResourceLocator(obj)), panel);
                    }
                } else {
                    messages.append("Unknown object type: " + transferable.toString());
                }
            } catch (Exception e) {
                log.error(e);
                if (e.getMessage() != null) {
                    messages.append(e.getMessage());
                }
            }

            if (messages != null && !messages.isEmpty()) {
                MessageUtils.showMessage(messages.getFormattedMessage());
            }

            igv.getMainFrame().repaint();
            event.dropComplete(true);
        }

        public boolean isDragAcceptable(DropTargetDragEvent event) { // usually, you
            // check the
            // available
            // data flavors
            // here
            // in this program, we accept all flavors
            return (event.getDropAction() & DnDConstants.ACTION_COPY_OR_MOVE) != 0;
        }

        public boolean isDropAcceptable(DropTargetDropEvent event) { // usually, you
            // check the
            // available
            // data flavors
            // here
            // in this program, we accept all flavors
            return (event.getDropAction() & DnDConstants.ACTION_COPY_OR_MOVE) != 0;
        }
    }
}
