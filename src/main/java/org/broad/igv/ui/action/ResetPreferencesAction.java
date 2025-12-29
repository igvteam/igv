/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * @author jrobinso
 */
public class ResetPreferencesAction extends MenuAction {

    // TODO -- The batch referenceFrame is likely to be used by many actions. Move this
    // member to a base class ?
    IGV mainFrame;

    public ResetPreferencesAction(String label, IGV mainFrame) {
        super(label);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.RESET_FACTORY_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                execute();

            }
        });
    }

    public void execute() {

        int status =
                JOptionPane.showConfirmDialog(
                        mainFrame.getMainFrame(),
                        "All application preferences will be reset to" +
                                " their default settings." +
                                "\n Continue?",
                        null,
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.PLAIN_MESSAGE,
                        null);

        if (status == JOptionPane.CANCEL_OPTION ||
                status == JOptionPane.CLOSED_OPTION) {
            return;
        }
        mainFrame.resetToFactorySettings();
    }
}
