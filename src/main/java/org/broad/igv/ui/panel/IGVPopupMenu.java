package org.broad.igv.ui.panel;

import org.broad.igv.ui.IGV;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.util.HashSet;

/**
 * Extension of JPopupMenu that provides an instance cache, and closeAll method.  Created so that popup
 * menus can be closed if the IGV main window looses focus.
 *
 * @author jrobinso
 * @date Jul 9, 2011
 */
public class IGVPopupMenu extends JPopupMenu {

    private static HashSet<IGVPopupMenu> instances = new HashSet();

    public IGVPopupMenu() {
        instances.add(this);
        addPopupMenuListener(new PopupMenuListener() {
            public void popupMenuWillBecomeVisible(PopupMenuEvent popupMenuEvent) {

            }

            public void popupMenuWillBecomeInvisible(PopupMenuEvent popupMenuEvent) {
                close();
            }

            public void popupMenuCanceled(PopupMenuEvent popupMenuEvent) {
                close();
            }

            private void close() {
                instances.remove(IGVPopupMenu.this);
            }

        });
    }

    public boolean includeStandardItems() {
        return true;  // Override to suppress save image, export, etc.
    }

    public static void closeAll() {
        synchronized (instances) {
            for (IGVPopupMenu inst : instances) {
                inst.setVisible(false);
            }
            instances.clear();
        }
    }


}
