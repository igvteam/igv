/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
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

    static Logger log = Logger.getLogger(NavigateRegionsMenuAction.class);
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
