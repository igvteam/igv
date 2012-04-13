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
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.SortDialog;
import org.broad.igv.ui.util.UIUtilities;

import java.awt.event.ActionEvent;
import java.util.List;

/**
 * @author jrobinso
 */
public class SortTracksMenuAction extends MenuAction {

    //static Logger log = Logger.getLogger(SortTracksMenuAction.class);
    IGV mainFrame;

    public SortTracksMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                doSortTrackByAttribute();
            }
        });

    }

    final public void doSortTrackByAttribute() {

        List<String> keys = AttributeManager.getInstance().getAttributeNames();
        Object availableSortKeys[] = keys.toArray();
        SortDialog dialog = new SortDialog(mainFrame.getMainFrame(), true, availableSortKeys);
        dialog.setVisible(true);

        if (dialog.isCanceled()) {
            return;
        }

        String[] selectedSortKeys = dialog.getSelectedSortKeys();
        if (selectedSortKeys != null) {
            IGV.getInstance().sortAllTracksByAttributes(selectedSortKeys, dialog.isAscending());
            mainFrame.getMainFrame().repaint();
        }


    }
}
