/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.ui.AttributeSelectionDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.UIUtilities;

import java.awt.event.ActionEvent;

/**
 * @author jrobinso
 */
public class GroupTracksMenuAction extends MenuAction {

    //static Logger log = LogManager.getLogger(GroupTracksMenuAction.class);
    IGV igv;

    public GroupTracksMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        doGroupBy();

    }

    final public void doGroupBy() {


        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {

                final AttributeSelectionDialog dlg = new AttributeSelectionDialog(
                        igv.getMainFrame(),
                        "Group");


                String currentSelection = IGV.getInstance().getGroupByAttribute();
                if (currentSelection == null) {
                    dlg.setSelectedIndex(0);
                } else {
                    dlg.setSelectedItem(currentSelection);
                }

                dlg.setVisible(true);

                if (!dlg.isCanceled()) {
                    String selectedAttribute = dlg.getSelected();
                    IGV.getInstance().setGroupByAttribute(selectedAttribute);
                    igv.repaint();

                }

            }

        });
    }
}
