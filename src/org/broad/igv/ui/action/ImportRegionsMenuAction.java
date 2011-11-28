/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;

import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

/**
 * @author jrobinso
 */
public class ImportRegionsMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(ImportRegionsMenuAction.class);
    IGV mainFrame;

    public ImportRegionsMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.IMPORT_REGION_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                File file = FileDialogUtils.chooseFile("Import regions of interest");
                if (file != null) {
                    readRegionsOfInterestFile(file);
                }
            }
        });
    }


    private void readRegionsOfInterestFile(File roiFile) {

        if (roiFile == null) {
            log.info("A blank Region of Interest import file was supplied!");
            return;
        }

        if (!roiFile.exists()) {
            MessageUtils.showMessage("Region of Interest export file not found!");
            return;
        }
        try {
            BufferedReader reader = null;
            int coordConvention = 0;

            try {
                reader = new BufferedReader(new FileReader(roiFile));
                while (true) {
                    String dataRecord = reader.readLine();
                    if (dataRecord == null) {
                        return;
                    } else if (dataRecord.startsWith("track")) {
                        // Skip track line
                        continue;
                    } else if (dataRecord.startsWith("#coords")) {
                        String[] tmp = dataRecord.split("=");
                        if (tmp.length > 1) {
                            try {
                                coordConvention = Integer.parseInt(tmp[1]);
                            } catch (NumberFormatException e) {
                                log.error("Error parsing coordinate convention direction for file: " + roiFile);
                            }
                        }
                    }
                    String[] data = dataRecord.split("\t");
                    if (data.length >= 3) {
                        try {
                            String name = data.length > 3 ? data[3] : null;
                            int start = Integer.parseInt(data[1]) - coordConvention;
                            int end = Integer.parseInt(data[2]);
                            RegionOfInterest regionOfInterest = new RegionOfInterest(data[0], start, end, name);
                            mainFrame.addRegionOfInterest(regionOfInterest);
                        } catch (NumberFormatException numberFormatException) {
                        }
                    }
                }
            } finally {

                if (reader != null) {
                    reader.close();
                }
                mainFrame.doRefresh();
            }
        } catch (Exception e) {
            log.error("Failed to write Region of Interest export file!", e);
        }
    }

}
