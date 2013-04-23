/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

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
                if (IGV.hasInstance()) {
                    IGV.getInstance().repaint();
                }
                instances.remove(IGVPopupMenu.this);
            }

        });
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
