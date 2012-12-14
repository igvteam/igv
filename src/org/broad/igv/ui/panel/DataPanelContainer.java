/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.MessageCollection;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 8, 2010
 */
public class DataPanelContainer extends TrackPanelComponent implements Paintable {

    private static Logger log = Logger.getLogger(DataPanelContainer.class);

    TrackPanel parent;

    public DataPanelContainer(TrackPanel trackPanel) {
        super(trackPanel);

        DropTarget target = new DropTarget(this, new FileDropTargetListener(trackPanel));
        setDropTarget(target);
        target.setActive(true);
        this.setLayout(new DataPanelLayout());
        this.parent = trackPanel;
        createDataPanels();

    }

    public void createDataPanels() {
        removeAll();
        for (ReferenceFrame f : FrameManager.getFrames()) {
            DataPanel dp = new DataPanel(f, this);
            add(dp);
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
     * @param g
     * @param rect
     */
    public void paintOffscreen(Graphics2D g, Rectangle rect) {

        // Get the components of the sort by X position.
        Component[] components = getComponents();
        Arrays.sort(components, new Comparator<Component>() {
            public int compare(Component component, Component component1) {
                return component.getX() - component1.getX();
            }
        });

        for (Component c : this.getComponents()) {

            if (c instanceof DataPanel) {
                Graphics2D g2d = (Graphics2D) g.create();
                g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                Rectangle clipRect = new Rectangle(c.getBounds());
                clipRect.height = rect.height;
                g2d.setClip(clipRect);
                g2d.translate(c.getX(), 0);
                ((DataPanel) c).paintOffscreen(g2d, rect);

            }
        }
        //super.paintBorder(g);
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

            MessageCollection messages = new MessageCollection();
            try {
                List<File> files = (List<File>) transferable.getTransferData(DataFlavor.javaFileListFlavor);

                for (File file : files) {
                    try {
                        ResourceLocator locator = new ResourceLocator(file.getAbsolutePath());
                        IGV.getInstance().load(locator, panel);
                    } catch (DataLoadException de) {
                        messages.append(de.getMessage());
                    }
                }
                String obj = transferable.getTransferData(DataFlavor.stringFlavor).toString();
                if (HttpUtils.isRemoteURL(obj)) {
                    IGV.getInstance().load(new ResourceLocator(obj), panel);
                }
                if (messages != null && !messages.isEmpty()) {
                    log.error(messages.getFormattedMessage());
                    MessageUtils.showMessage(messages.getFormattedMessage());
                }
            } catch (Exception e) {
                String obj = null;
                try {
                    obj = transferable.getTransferData(DataFlavor.stringFlavor).toString();
                    if (HttpUtils.isRemoteURL(obj)) {
                        IGV.getInstance().load(new ResourceLocator(obj), panel);
                    }
                } catch (Exception e1) {
                    log.error(e1);
                    if (messages != null && !messages.isEmpty()) {
                        MessageUtils.showMessage(messages.getFormattedMessage());
                    }
                }
            }
            IGV.getMainFrame().repaint();
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
