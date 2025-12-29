package org.broad.igv.ui.panel;

import javax.swing.*;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceMotionListener;

/**
 * <p>Used by both the draggable class and the target for negotiating data.</p>
 * <p>Note that this should be set for both the draggable object and the drop target.</p>
 *
 * @author besmit
 */
class DragAndDropTransferHandler extends TransferHandler implements DragSourceMotionListener {

    public DragAndDropTransferHandler() {
        super();
    }

    /**
     * <p>This creates the Transferable object. In our case, RandomDragAndDropPanel implements Transferable, so this requires only a type cast.</p>
     *
     * @param c
     * @return
     */
    @Override()
    public Transferable createTransferable(JComponent c) {

        // TaskInstancePanel implements Transferable
        if (c instanceof HeaderPanel) {
            Transferable tip = (HeaderPanel) c;
            return tip;
        }

        // Not found
        return null;
    }

    public void dragMouseMoved(DragSourceDragEvent dsde) {
    }

    /**
     * <p>This is queried to see whether the component can be copied, moved, both or neither. We are only concerned with copying.</p>
     *
     * @param c
     * @return
     */
    @Override()
    public int getSourceActions(JComponent c) {

        if (c instanceof HeaderPanel) {
            return TransferHandler.MOVE;
        }

        return TransferHandler.NONE;
    }
} // DragAndDropTransferHandler
