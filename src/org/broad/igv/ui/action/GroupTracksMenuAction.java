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

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                doGroupBy();
            }
        });

    }

    final public void doGroupBy() {

        final AttributeSelectionDialog dlg = new AttributeSelectionDialog(mainFrame.getMainFrame(), true);

        List<String> attributeKeys = AttributeManager.getInstance().getAttributeNames();


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
}
