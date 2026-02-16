package org.igv.ui.panel;

import org.igv.ui.IGV;

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
                dtde.rejectDrop();
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

        // Determine if drop is in top or bottom half of target panel
        int dropYLoc = dtde.getLocation().y;
        boolean before = dropYLoc < targetPanel.getHeight() / 2;

        // Reorder the panels
        MainPanel mainPanel = IGV.getInstance().getMainPanel();
        List<TrackPanel> panels = mainPanel.getTrackPanels();
        List<String> orderedNames = new ArrayList<>(panels.size());

        // Remove the dropped panel from its current position
        for (TrackPanel panel : panels) {
            if (panel != droppedPanel) {
                if (panel == targetPanel) {
                    if (before) {
                        orderedNames.add(droppedPanel.getName());
                        orderedNames.add(panel.getName());
                    } else {
                        orderedNames.add(panel.getName());
                        orderedNames.add(droppedPanel.getName());
                    }
                } else {
                    orderedNames.add(panel.getName());
                }
            }
        }

        // Reorder the panels in the MainPanel
        mainPanel.reorderPanels(orderedNames);

        dtde.dropComplete(true);
    }
}

