package org.broad.igv.ui.util;

import org.broad.igv.ui.UIConstants;

import javax.swing.event.MouseInputAdapter;
import java.awt.event.MouseEvent;

public class IGVMouseInputAdapter extends MouseInputAdapter {

    MouseEvent mouseDown;
    long lastClickTime = 0;
    int clickCount = 1;

    @Override
    public void mousePressed(final MouseEvent e) {
        if (e.isPopupTrigger()) {
            // ignore
        } else {
            mouseDown = e;
        }
    }


    @Override
    public void mouseReleased(MouseEvent e) {
        if (e.isPopupTrigger()) {
            // ignore
        } else {
            if (mouseDown != null && distance(mouseDown, e) < 5) {
                long time = System.currentTimeMillis();
                if (time - lastClickTime < UIConstants.getDoubleClickInterval()) {
                    clickCount++;
                } else {
                    clickCount = 1;
                }
                lastClickTime = time;
                igvMouseClicked(new IGVMouseEvent(e, clickCount));
            }
        }
    }

    public void igvMouseClicked(MouseEvent e) {
    }

    private double distance(MouseEvent e1, MouseEvent e2) {
        double dx = e1.getX() - e2.getX();
        double dy = e1.getY() - e2.getY();
        return Math.sqrt(dx * dx + dy * dy);
    }
}

class IGVMouseEvent extends MouseEvent {

    IGVMouseEvent(MouseEvent e, int clickCount) {
        super(e.getComponent(), e.getID(), e.getWhen(), e.getModifiersEx(), e.getX(), e.getY(), clickCount, e.isPopupTrigger(), e.getButton());
    }

}