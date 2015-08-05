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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.AttributeSelectionDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class GroupTracksMenuAction extends MenuAction {

    //static Logger log = Logger.getLogger(GroupTracksMenuAction.class);
    IGV mainFrame;

    public GroupTracksMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        doGroupBy();

    }

    final public void doGroupBy() {


        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {

                final AttributeSelectionDialog dlg = new AttributeSelectionDialog(mainFrame.getMainFrame(), true);

                List<String> attributeKeys = AttributeManager.getInstance().getVisibleAttributes();


                // Sorting disabled -- order will match the order in the panel.  If sorting is desired make a copy
                // of the array so the panel is not affected.

                //if (attributeKeys != null) {
                //    Collections.sort(attributeKeys,
                //            AttributeManager.getInstance().getAttributeComparator());
                //}

                ArrayList<String> selections = new ArrayList(attributeKeys);

                selections.add(0, "None");
                String[] selArray = selections.toArray(new String[]{});

                dlg.setModel(new javax.swing.DefaultComboBoxModel(selArray));
                dlg.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

                String currentSelection = IGV.getInstance().getGroupByAttribute();
                if (currentSelection == null) {
                    dlg.setSelectedIndex(0);
                } else {
                    dlg.setSelectedItem(currentSelection);
                }

                dlg.setVisible(true);

                if (!dlg.isCanceled()) {
                    int selIndex = dlg.getSelectedIndex();
                    String selectedAttribute = (selIndex == 0 ? null : selArray[selIndex]);
                    IGV.getInstance().setGroupByAttribute(selectedAttribute);
                    mainFrame.doRefresh();

                }

            }

        });
    }
}
