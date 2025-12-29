/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.action;

import org.igv.track.AttributeManager;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.util.SortDialog;
import org.igv.ui.util.UIUtilities;

import java.awt.event.ActionEvent;
import java.util.List;

/**
 * @author jrobinso
 */
public class SortTracksMenuAction extends MenuAction {

    //static Logger log = LogManager.getLogger(SortTracksMenuAction.class);
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
