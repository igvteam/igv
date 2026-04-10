package org.igv.ui.panel;

import javax.swing.*;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceMotionListener;

/**
 * Handler for moving TrackPanels via drag and drop.
 * Used by both the draggable TrackPanel and the drop target for negotiating data.
 */
class TrackPanelTransferHandler extends TransferHandler implements DragSourceMotionListener {

    public TrackPanelTransferHandler() {
        super();
    }

    /**
     * Creates the Transferable object. TrackPanel implements Transferable,
     * so this requires only a type cast.
     */
    @Override
    public Transferable createTransferable(JComponent c) {
        if (c instanceof TrackPanel) {
            return (TrackPanel) c;
        }
        return null;
    }

    @Override
    public void dragMouseMoved(DragSourceDragEvent dsde) {
    }

    /**
     * Returns the supported source actions. We only support MOVE.
     */
    @Override
    public int getSourceActions(JComponent c) {
        if (c instanceof TrackPanel) {
            return TransferHandler.MOVE;
        }
        return TransferHandler.NONE;
    }
}

