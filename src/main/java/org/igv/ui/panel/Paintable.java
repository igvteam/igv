package org.igv.ui.panel;

import javax.swing.*;
import java.awt.*;

/**
 * @author jrobinso
 * @date Dec 8, 2010
 */
public interface Paintable  {

    void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch);

    int getSnapshotHeight(boolean batch);

}
