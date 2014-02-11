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

    int height;
    int trackHeight = 60;
    int vmargin = 0;

    @Override
    public void doLayout() {
        height = 0;
        synchronized (this.getTreeLock()) {
            final int width = getWidth();
            for (Component c : this.getComponents()) {
                c.setBounds(0, height, width, trackHeight - 2*vmargin);
                height += trackHeight + vmargin;
            }
        }
    }

    @Override
    public int getHeight() {
        return height;
    }

    @Override
    public Component add(Component comp) {
        int sz = (this.getComponents().length + 1) * trackHeight;
        setPreferredSize(new Dimension(getWidth(), sz));
        return super.add(comp);
    }

    public void setVmargin(int vmargin) {
        this.vmargin = vmargin;
    }

}
