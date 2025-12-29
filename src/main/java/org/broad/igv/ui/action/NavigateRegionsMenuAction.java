/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.logging.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.RegionNavigatorDialog;
import org.broad.igv.ui.util.UIUtilities;

import java.awt.event.ActionEvent;

/**
 * @author Damon May
 *         Menu action that brings up a region navigation dialog box
 */
public class NavigateRegionsMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(NavigateRegionsMenuAction.class);
    IGV mainFrame;


    public NavigateRegionsMenuAction(String label, IGV mainFrame) {
        super(label, null);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.REGION_NAVIGATOR_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                RegionNavigatorDialog regionNavDialog = RegionNavigatorDialog.getOrCreateInstance(mainFrame.getMainFrame());
                regionNavDialog.setVisible(true);
            }
        });
    }
}
