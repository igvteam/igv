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
 * TrackPanel.java
 *
 * Created on Sep 5, 2007, 4:09:39 PM
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.ui.panel;


import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.SwitchingLabelUI;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;


/**
 * The drag & drop code was modified from the excellent example of Bryan E. Smith.
 *
 * @author jrobinso
 */
public class HeaderPanel extends JPanel implements Transferable {

    ReferenceFrame frame;
    private JLabel label;
    public static final Color BUTTON_BACKGROUND = new Color(230, 240, 250);
    static DataFlavor dragAndDropPanelDataFlavor;
    private CytobandPanel cytobandPanel;
    private RulerPanel rulerPanel;
    private RegionOfInterestPanel regionOfInterestPane;
    private JPanel geneListPanel;


    public HeaderPanel(ReferenceFrame frame) {
        this.frame = frame;
        init();
    }

    private void init() {


        setBackground(new java.awt.Color(255, 255, 255));
        setMinimumSize(new java.awt.Dimension(700, 0));
        setPreferredSize(new java.awt.Dimension(0, 0));
        setLayout(new java.awt.BorderLayout());

        if (FrameManager.isGeneListMode()) {

            geneListPanel = new JPanel();
            geneListPanel.setMinimumSize(new java.awt.Dimension(700, 0));
            geneListPanel.setPreferredSize(new java.awt.Dimension(0, 0));
            geneListPanel.setLayout(new java.awt.BorderLayout());


            label = new JLabel(frame.getName());
            label.setForeground(Color.blue);
            label.setUI(new SwitchingLabelUI(10));
            label.setToolTipText(frame.getName());
            label.setPreferredSize(new Dimension(500, 80));

            //label.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            //label.addMouseListener(mouseAdapter);


            setSize(400, 100);
            setVisible(true);

            cytobandPanel = new CytobandPanel(frame, false);
            cytobandPanel.setBackground(new java.awt.Color(255, 255, 255));
            cytobandPanel.setPreferredSize(new java.awt.Dimension(0, 20));
            cytobandPanel.setRequestFocusEnabled(false);
            cytobandPanel.setLayout(null);

            final MouseAdapter mouseAdapter = new MouseAdapter() {

                boolean isDragging = false;
                Point mousePressPoint;

                @Override
                public void mouseClicked(MouseEvent mouseEvent) {
                    if (mouseEvent.getClickCount() > 1) {
                        IGV.getInstance().setDefaultFrame(frame.getName());
                    }
                }


                public void mousePressed(MouseEvent evt) {
                    if (evt.isPopupTrigger()) {
                        getPopupMenu(HeaderPanel.this, frame).show(HeaderPanel.this, evt.getX(), evt.getY());
                    }
                }

                @Override
                public void mouseReleased(MouseEvent evt) {

                    if (evt.isPopupTrigger()) {
                        getPopupMenu(HeaderPanel.this, frame).show(HeaderPanel.this, evt.getX(), evt.getY());
                    }
                    isDragging = false;

                }


                @Override()
                public void mouseDragged(MouseEvent e) {

                    if (isDragging) {
                        return;
                    }

                    isDragging = true;
                    JComponent c = HeaderPanel.this;
                    TransferHandler handler = c.getTransferHandler();
                    if (handler != null) {
                        handler.exportAsDrag(c, e, TransferHandler.MOVE);
                    }

                }
            };

            cytobandPanel.addMouseListener(mouseAdapter);
            cytobandPanel.addMouseMotionListener(mouseAdapter);

            label.addMouseListener(mouseAdapter);
            label.addMouseMotionListener(mouseAdapter);

            this.addMouseListener(mouseAdapter);
            this.addMouseMotionListener(mouseAdapter);

            geneListPanel.add(cytobandPanel, BorderLayout.CENTER);
            geneListPanel.add(label, BorderLayout.SOUTH);
            add(geneListPanel);

            this.setTransferHandler(new DragAndDropTransferHandler());
            // Create the listener to do the work when dropping on this object!
            this.setDropTarget(new DropTarget(this, new HeaderDropTargetListener(this)));


        } else {

            JPanel panel = new JPanel();
            setBorder(javax.swing.BorderFactory.createLineBorder(Color.gray));
            panel.setBackground(new java.awt.Color(255, 255, 255));
            panel.setMinimumSize(new java.awt.Dimension(700, 0));
            panel.setPreferredSize(new java.awt.Dimension(0, 0));
            panel.setLayout(new java.awt.BorderLayout());

            cytobandPanel = new CytobandPanel(frame);
            cytobandPanel.setBackground(new java.awt.Color(255, 255, 255));
            cytobandPanel.setPreferredSize(new java.awt.Dimension(0, 50));
            cytobandPanel.setRequestFocusEnabled(false);
            cytobandPanel.setLayout(null);
            panel.add(cytobandPanel, java.awt.BorderLayout.NORTH);

            rulerPanel = new RulerPanel(frame);
            rulerPanel.setBackground(new java.awt.Color(255, 255, 255));
            rulerPanel.setLayout(null);
            panel.add(rulerPanel, java.awt.BorderLayout.CENTER);

            regionOfInterestPane = new RegionOfInterestPanel(frame);
            regionOfInterestPane.setBackground(new java.awt.Color(255, 255, 255));
            regionOfInterestPane.setMinimumSize(new java.awt.Dimension(0, 13));


            panel.add(regionOfInterestPane, java.awt.BorderLayout.SOUTH);

            add(panel);
        }


    }


    // TODO -- this is a partial copy of the RegionOfInterestPanel method.  Refactor to share


    protected JPopupMenu getPopupMenu(final HeaderPanel parent, final ReferenceFrame frame) {

        int start = (int) frame.getOrigin();
        int end = (int) frame.getEnd();
        final RegionOfInterest roi = new RegionOfInterest(frame.getChrName(), start, end, "");

        JPopupMenu popupMenu = new RegionMenu(roi, frame, "Panel: " + frame.getName());


        // Zoom items
        popupMenu.addSeparator();
        TrackMenuUtils.addZoomItems(popupMenu, frame);
        popupMenu.addSeparator();
        JMenuItem item1 = new JMenuItem("Switch to standard view");
        item1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                IGV.getInstance().setDefaultFrame(frame.getName());

            }
        });
        popupMenu.add(item1);

        popupMenu.addSeparator();
        JMenuItem item = new JMenuItem("Remove panel");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                parent.removeFrame(frame);
            }
        });
        popupMenu.add(item);


        return popupMenu;
    }

    private void removeFrame(ReferenceFrame frame) {
        FrameManager.removeFrame(frame);
        java.util.List<ReferenceFrame> remainingFrames = FrameManager.getFrames();
        if (remainingFrames.size() == 1) {
            IGV.getInstance().setDefaultFrame(remainingFrames.get(0).getName());
        } else {
            IGV.getInstance().resetFrames();
        }
    }

    @Override

    public void setBackground(Color color) {
        super.setBackground(color);
        if (cytobandPanel != null) cytobandPanel.setBackground(color);
        if (rulerPanel != null) rulerPanel.setBackground(color);
        if (regionOfInterestPane != null) regionOfInterestPane.setBackground(color);
        if (geneListPanel != null) geneListPanel.setBackground(color);
    }


    /**
     * <p>Returns (creating, if necessary) the DataFlavor representing RandomDragAndDropPanel</p>
     *
     * @return
     */
    public static DataFlavor getDragAndDropPanelDataFlavor() throws Exception {
        // Lazy load/create the flavor
        if (dragAndDropPanelDataFlavor == null) {
            dragAndDropPanelDataFlavor = new DataFlavor(DataFlavor.javaJVMLocalObjectMimeType +
                    ";class=org.broad.igv.ui.panel.HeaderPanel");
        }

        return dragAndDropPanelDataFlavor;
    }

    //private static final Cursor droppableCursor = Cursor.getPredefinedCursor(Cursor.HAND_CURSOR);
    //private static final Cursor notDroppableCursor = Cursor.getDefaultCursor();

    /**
     * <p>Listens for drops and performs the updates.</p>
     * <p>The real magic behind the drop!</p>
     */
    class HeaderDropTargetListener implements DropTargetListener {

        private final HeaderPanel rootPanel;

        /**
         * <p>Two cursors with which we are primarily interested while dragging:</p>
         * <ul>
         * <li>Cursor for droppable condition</li>
         * <li>Cursor for non-droppable consition</li>
         * </ul>
         * <p>After drop, we manually change the cursor back to default, though does this anyhow -- just to be complete.</p>
         */

        public HeaderDropTargetListener(HeaderPanel sheet) {
            this.rootPanel = sheet;
        }

        // Could easily find uses for these, like cursor changes, etc.

        public void dragEnter(DropTargetDragEvent dtde) {
        }

        public void dragOver(DropTargetDragEvent dtde) {
            // if (!this.rootPanel.getCursor().equals(droppableCursor)) {
            //     this.rootPanel.setCursor(droppableCursor);
            // }
        }

        public void dropActionChanged(DropTargetDragEvent dtde) {
        }

        public void dragExit(DropTargetEvent dte) {
            // this.rootPanel.setCursor(notDroppableCursor);
        }

        /**
         * <p>The user drops the item. Performs the drag and drop calculations and layout.</p>
         *
         * @param dtde
         */
        public void drop(DropTargetDropEvent dtde) {

            // Done with cursors, dropping
            //this.rootPanel.setCursor(Cursor.getDefaultCursor());

            // Just going to grab the expected DataFlavor to make sure
            // we know what is being dropped
            DataFlavor dragAndDropPanelFlavor = null;
            Object transferableObj = null;

            try {
                // Grab expected flavor
                dragAndDropPanelFlavor = HeaderPanel.getDragAndDropPanelDataFlavor();

                Transferable transferable = dtde.getTransferable();

                // What does the Transferable support
                if (transferable.isDataFlavorSupported(dragAndDropPanelFlavor)) {
                    transferableObj = dtde.getTransferable().getTransferData(dragAndDropPanelFlavor);
                }

            } catch (Exception ex) { /* nope, not the place */ }

            // If didn't find an item, bail
            if (transferableObj == null) {
                return;
            }

            // Cast it to the panel. By this point, we have verified it is a HeaderPanel
            HeaderPanel droppedPanel = (HeaderPanel) transferableObj;
            ReferenceFrame droppedFrame = droppedPanel.frame;
            if (droppedFrame == frame) {
                IGV.getInstance().resetFrames();
            } else {
                final int dropXLoc = dtde.getLocation().x;
                boolean before = dropXLoc < getWidth() / 2;


                // Find the index for the drop
                java.util.List<ReferenceFrame> panels = FrameManager.getFrames();
                java.util.List<ReferenceFrame> orderedPanels = new ArrayList(panels.size());
                panels.remove(droppedFrame);

                boolean dropAdded = false;


                for (ReferenceFrame frame : panels) {
                    if (HeaderPanel.this.frame == frame) {
                        if (before) {
                            orderedPanels.add(droppedFrame);
                            orderedPanels.add(frame);
                        } else {
                            orderedPanels.add(frame);
                            orderedPanels.add(droppedFrame);
                        }
                        dropAdded = true;
                    } else {
                        orderedPanels.add(frame);
                    }
                }


                // Request relayout contents, or else won't update GUI following drop.
                // Will add back in the order to which we just sorted
                FrameManager.setFrames(orderedPanels);
                IGV.getInstance().resetFrames();
            }
        }
    } // HeaderDropTargetListener

    /**
     * <p>One of three methods defined by the Transferable interface.</p>
     * <p>If multiple DataFlavor's are supported, can choose what Object to return.</p>
     * <p>In this case, we only support one: the actual JPanel.</p>
     * <p>Note we could easily support more than one. For example, if supports text and drops to a JTextField, could return the label's text or any arbitrary text.</p>
     *
     * @param flavor
     * @return
     */
    public Object getTransferData(DataFlavor flavor) {

        DataFlavor thisFlavor = null;

        try {
            thisFlavor = getDragAndDropPanelDataFlavor();
        } catch (Exception ex) {
            System.err.println("Problem lazy loading: " + ex.getMessage());
            ex.printStackTrace(System.err);
            return null;
        }

        // For now, assume wants this class... see loadDnD
        if (thisFlavor != null && flavor.equals(thisFlavor)) {
            return this;
        }

        return null;
    }

    /**
     * <p>One of three methods defined by the Transferable interface.</p>
     * <p>Returns supported DataFlavor. Again, we're only supporting this actual Object within the JVM.</p>
     * <p>For more information, see the JavaDoc for DataFlavor.</p>
     *
     * @return
     */
    public DataFlavor[] getTransferDataFlavors() {

        DataFlavor[] flavors = {null};
        try {
            flavors[0] = getDragAndDropPanelDataFlavor();
        } catch (Exception ex) {
            System.err.println("Problem lazy loading: " + ex.getMessage());
            ex.printStackTrace(System.err);
            return null;
        }

        return flavors;
    }

    /**
     * <p>One of three methods defined by the Transferable interface.</p>
     * <p>Determines whether this object supports the DataFlavor. In this case, only one is supported: for this object itself.</p>
     *
     * @param flavor
     * @return True if DataFlavor is supported, otherwise false.
     */
    public boolean isDataFlavorSupported(DataFlavor flavor) {

        DataFlavor[] flavors = {null};
        try {
            flavors[0] = getDragAndDropPanelDataFlavor();
        } catch (Exception ex) {
            System.err.println("Problem lazy loading: " + ex.getMessage());
            ex.printStackTrace(System.err);
            return false;
        }

        for (DataFlavor f : flavors) {
            if (f.equals(flavor)) {
                return true;
            }
        }

        return false;
    }

}

