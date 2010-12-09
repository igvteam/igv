/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.RegionNavigatorDialog;
import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.Collection;

/**
 * @author Damon May
 * Menu action that brings up a region navigation dialog box
 */
public class NavigateRegionsMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(NavigateRegionsMenuAction.class);
    IGVMainFrame mainFrame;



    public NavigateRegionsMenuAction(String label, IGVMainFrame mainFrame) {
        super(label, null);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.NAVIGATE_REGION_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
//                Collection<RegionOfInterest> regions = IGVMainFrame.getInstance().getSession().getAllRegionsOfInterest();
//                if (regions == null || regions.isEmpty())
//                {
//                    //todo dhmay -- I don't fully understand this call.  Clean this up.
//                    JOptionPane.showMessageDialog(mainFrame, "No regions of interest are currently loaded",
//                            "Error", JOptionPane.INFORMATION_MESSAGE);
//                }
//                else
//                {
                    if (RegionNavigatorDialog.getActiveInstance() == null)
                    {
                        RegionNavigatorDialog regionNavDialog = new RegionNavigatorDialog(mainFrame);
                    }
                    RegionNavigatorDialog.getActiveInstance().setVisible(true);
//                }
            }
        });
    }
}
