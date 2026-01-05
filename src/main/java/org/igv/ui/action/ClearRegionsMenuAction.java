package org.igv.ui.action;

import org.igv.logging.*;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.util.UIUtilities;

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
