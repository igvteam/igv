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
