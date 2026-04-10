package org.igv.ui.panel;

import org.igv.Globals;

import javax.swing.*;
import java.awt.*;

/**
 * Name header panel.  Displayed in upper left corner
 *
 * @author jrobinso
 */
public class NameHeaderPanel extends JPanel implements Paintable {


    private final boolean darkMode;

    public NameHeaderPanel() {
        this.darkMode = Globals.isDarkMode();
        setBorder(null);
        //if(darkMode){
        //    setBackground(UIManager.getColor("Panel.background"));
        //}
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {
        paintComponent(g);
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
    }
}
