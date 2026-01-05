package org.igv.ui.action;

import org.igv.logging.*;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.panel.RegionNavigatorDialog;
import org.igv.ui.util.UIUtilities;

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
