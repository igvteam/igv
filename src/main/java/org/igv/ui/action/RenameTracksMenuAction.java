/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.action;

import org.igv.prefs.Constants;
import org.igv.session.Session;
import org.igv.ui.AttributeSelectionDialog;
import org.igv.ui.IGV;
import org.igv.ui.util.UIUtilities;

import java.awt.event.ActionEvent;

/**
 * @author jrobinso
 */
public class RenameTracksMenuAction extends MenuAction {

    IGV igv;

    public RenameTracksMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        UIUtilities.invokeOnEventThread(() -> {

            final AttributeSelectionDialog dlg = new AttributeSelectionDialog(igv.getMainFrame(), "Rename");
            dlg.setVisible(true);

            if (!dlg.isCanceled()) {
                Session currentSession = IGV.getInstance().getSession();
                String selectedAttribute = dlg.getSelected();
                if (selectedAttribute == null) {
                    currentSession.removePreference(Constants.TRACK_ATTRIBUTE_NAME_KEY);
                } else {
                    currentSession.setPreference(Constants.TRACK_ATTRIBUTE_NAME_KEY, selectedAttribute);
                }
                igv.repaint();
            }
        });
    }

}
