/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
