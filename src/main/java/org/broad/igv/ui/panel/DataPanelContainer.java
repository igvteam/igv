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

package org.broad.igv.ui.panel;

import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.load.HubGenomeLoader;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.MessageCollection;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
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
            try {
                if (transferable.isDataFlavorSupported(DataFlavor.javaFileListFlavor)) {
                    List<File> files = (List<File>) transferable.getTransferData(DataFlavor.javaFileListFlavor);
                    if (files != null && files.size() > 0) {
                        List<ResourceLocator> locators = ResourceLocator.getLocators(files);
                        for (ResourceLocator locator : locators) {
                            try {
                                IGV.getInstance().load(locator, panel);
                            } catch (DataLoadException de) {
                                messages.append(de.getMessage());
                            }
                        }
                    }
                } else if (transferable.isDataFlavorSupported(DataFlavor.stringFlavor)) {
                    String obj = transferable.getTransferData(DataFlavor.stringFlavor).toString();

                    // Check for genomes, sessions, etc firs
                    if(HubGenomeLoader.isHubURL(obj)) {
                       LongRunningTask.submit(() -> {
                           try {
                               GenomeManager.getInstance().loadGenome(obj);
                           } catch (IOException e) {
                               MessageUtils.showMessage("Error loading track hub: " + e.getMessage());
                           }
                       });
                    }
                    else {
                        IGV.getInstance().load(new ResourceLocator(obj), panel);
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

            IGV.getInstance().getMainFrame().repaint();
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
