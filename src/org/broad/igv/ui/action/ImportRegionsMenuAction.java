/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
            log.info("No regions found in file");
            return;
        }

        if (!roiFile.exists()) {
            MessageUtils.showMessage(roiFile.getAbsolutePath() + "  not found");
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
                    if (data.length >= 2) {
                        try {
                            String name = data.length > 3 ? data[3] : null;
                            int start = Integer.parseInt(data[1]) - coordConvention;
                            int end = data.length > 2 ? Integer.parseInt(data[2]) : start + 1;
                            RegionOfInterest regionOfInterest = new RegionOfInterest(data[0], start, end, name);
                            mainFrame.addRegionOfInterest(regionOfInterest);
                        } catch (NumberFormatException numberFormatException) {
                            log.error("Error importing regions of interest", numberFormatException);
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
