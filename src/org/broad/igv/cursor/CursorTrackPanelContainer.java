package org.broad.igv.cursor;

import javax.swing.*;
import java.awt.*;
import java.io.Serializable;

/**
 * @author jrobinso
 *         Date: 1/26/14
 *         Time: 5:11 PM
 */
public class CursorTrackPanelContainer extends JPanel implements Serializable {

    int trackHeight = 40;


    @Override
    public void doLayout() {
        synchronized (this.getTreeLock()) {
            int top = 0;
            final int width = getWidth();
            for (Component c : this.getComponents()) {
                c.setBounds(0, top, width, trackHeight);
                top += trackHeight;

            }
        }
    }

//    @Override
//    public int getHeight() {
//        return this.getComponents().length * trackHeight;
//    }


}
