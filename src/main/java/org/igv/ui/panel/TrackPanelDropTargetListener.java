package org.igv.ui.panel;

import org.igv.ui.IGV;

import javax.swing.SwingUtilities;
import java.awt.Component;
import java.awt.Point;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Listens for drops on TrackPanels and performs the reordering.
 */
class TrackPanelDropTargetListener implements DropTargetListener {

    private final TrackPanel targetPanel;

    public TrackPanelDropTargetListener(TrackPanel panel) {
        this.targetPanel = panel;
    }

    @Override
    public void dragEnter(DropTargetDragEvent dtde) {
    }

    @Override
    public void dragOver(DropTargetDragEvent dtde) {
    }

    @Override
    public void dropActionChanged(DropTargetDragEvent dtde) {
    }

    @Override
    public void dragExit(DropTargetEvent dte) {
    }

    /**
     * The user drops the item. Performs the drag and drop calculations and layout.
     */
    @Override
    public void drop(DropTargetDropEvent dtde) {
        DataFlavor trackPanelFlavor = null;
        Object transferableObj = null;

        try {
            trackPanelFlavor = TrackPanel.getTrackPanelDataFlavor();
            Transferable transferable = dtde.getTransferable();

            if (transferable.isDataFlavorSupported(trackPanelFlavor)) {
                dtde.acceptDrop(DnDConstants.ACTION_MOVE);
                transferableObj = transferable.getTransferData(trackPanelFlavor);
            } else {
                // Not a TrackPanel drag — delegate to MainPanel for file/URL drops
                MainPanel mainPanel = IGV.getInstance().getMainPanel();
                mainPanel.drop(dtde);
                return;
            }
        } catch (Exception ex) {
            dtde.rejectDrop();
            return;
        }

        if (transferableObj == null) {
            dtde.dropComplete(false);
            return;
        }

        TrackPanel droppedPanel = (TrackPanel) transferableObj;

        // If dropped on itself, do nothing
        if (droppedPanel == targetPanel) {
            dtde.dropComplete(true);
            return;
        }

        // Determine if drop is in top or bottom half of the *visible* track area.
        // The drop may have landed on any descendant of TrackPanelScrollPane (track
        // panel, name panel, data panel, drag handle, …) so convert the drop point
        // into the scroll pane's coordinate space and compare against its height —
        // that is always the visible track region regardless of scroll offset.
        Component dropComp = dtde.getDropTargetContext().getComponent();
        Point dropPoint = dtde.getLocation();
        TrackPanelScrollPane scrollPane = targetPanel.getScrollPane();
        boolean before;
        if (scrollPane != null && dropComp != null) {
            Point inScrollPane = SwingUtilities.convertPoint(dropComp, dropPoint, scrollPane);
            before = inScrollPane.y < scrollPane.getHeight() / 2;
        } else {
            before = dropPoint.y < targetPanel.getHeight() / 2;
        }

        // Reorder the panels using direct scroll pane references
        MainPanel mainPanel = IGV.getInstance().getMainPanel();
        List<TrackPanel> panels = mainPanel.getTrackPanels();
        List<TrackPanelScrollPane> orderedPanes = new ArrayList<>(panels.size());

        // Remove the dropped panel from its current position and insert at target
        for (TrackPanel panel : panels) {
            if (panel == droppedPanel) continue;
            if (panel == targetPanel) {
                if (before) {
                    orderedPanes.add(droppedPanel.getScrollPane());
                    orderedPanes.add(panel.getScrollPane());
                } else {
                    orderedPanes.add(panel.getScrollPane());
                    orderedPanes.add(droppedPanel.getScrollPane());
                }
            } else {
                orderedPanes.add(panel.getScrollPane());
            }
        }

        // Reorder the panels in the MainPanel
        mainPanel.reorderPanels(orderedPanes);

        // Update the moved track's order property to reflect its new position
        mainPanel.updateMovedTrackOrder(droppedPanel);

        dtde.dropComplete(true);
    }
}

