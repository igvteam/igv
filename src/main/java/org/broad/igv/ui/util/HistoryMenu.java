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
