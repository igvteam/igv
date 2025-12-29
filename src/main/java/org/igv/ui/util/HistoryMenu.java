package org.igv.ui.util;

import org.igv.ui.IGV;
import org.igv.session.History;
import org.igv.ui.MenuSelectedListener;

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


        this.addMenuListener(new MenuSelectedListener() {
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
        });
    }


}
