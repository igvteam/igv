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

package org.broad.igv.ui.util;

import org.broad.igv.ui.IGV;
import org.broad.igv.session.History;

import javax.swing.*;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 1, 2010
 */
public class HistoryMenu extends JMenu {
    static int MAX_ITEMS = 30;
    private JMenuItem backItem;
    private JMenuItem forwardItem;
    private JMenuItem clearAllItem;

    public HistoryMenu() {
        this("History");
    }


    public HistoryMenu(String name) {
        super(name);
        final History history = IGV.getInstance().getSession().getHistory();

        clearAllItem = new JMenuItem("Clear all");
        clearAllItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                history.clear();
            }
        });

        //backItem = new JMenuItem("<html>Back&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>Alt &rarr;");
        backItem = new JMenuItem("Back          Alt+Arrow");
        backItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                history.back();
            }
        });

        forwardItem = new JMenuItem("Forward     Alt+Arrow");
        forwardItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                history.forward();
            }
        });


        this.addMenuListener(new MenuListener() {
            public void menuSelected(MenuEvent menuEvent) {

                final History history = IGV.getInstance().getSession().getHistory();


                List<History.Entry> allLoci = IGV.getInstance().getSession().getAllHistory();
                
                boolean hasBack = history.peekBack() != null;
                boolean hasForward = history.peekForward() != null;
                backItem.setEnabled(hasBack);
                forwardItem.setEnabled(hasForward);
                clearAllItem.setEnabled(allLoci.size() > 0);

                // Update history list
                removeAll();

                add(backItem);
                add(forwardItem);
                addSeparator();

                int nItems = 0;
                // Do in reverse order

                for (int idx = allLoci.size() - 1; idx >= 0; idx--) {
                    final History.Entry s = allLoci.get(idx);
                    final JMenuItem item = new JMenuItem(s.toString());
                    item.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent actionEvent) {
                            IGV.getInstance().getSession().getHistory().processItem(s);
                        }
                    });
                    add(item);
                    if (nItems++ > MAX_ITEMS) {
                        break;
                    }
                }


                addSeparator();
                add(clearAllItem);


            }

            public void menuDeselected(MenuEvent menuEvent) {
                //To change body of implemented methods use File | Settings | File Templates.
            }

            public void menuCanceled(MenuEvent menuEvent) {
                //To change body of implemented methods use File | Settings | File Templates.
            }
        });
    }


}
