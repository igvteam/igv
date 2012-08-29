/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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

    static HashSet<IGVPopupMenu> instances = new HashSet();


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
                    IGV.getInstance().clearSelections();
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
