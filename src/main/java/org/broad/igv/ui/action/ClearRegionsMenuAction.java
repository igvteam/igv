/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.logging.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * @author jrobinso
 */
public class ClearRegionsMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(ClearRegionsMenuAction.class);
    IGV igv;

    public ClearRegionsMenuAction(String label, IGV igv) {
        super(label, null);
        this.igv = igv;
        setToolTipText(UIConstants.EXPORT_REGION_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                int choice = JOptionPane.showConfirmDialog(
                        igv.getMainFrame(),
                        "This action will clear all regions of interest.  Continue?",
                        "Clear Regions",
                        JOptionPane.YES_NO_OPTION);
                if (choice == JOptionPane.YES_OPTION) {
                    igv.getSession().clearRegionsOfInterest();
                    igv.repaint();
                }
            }
        });
    }
}
