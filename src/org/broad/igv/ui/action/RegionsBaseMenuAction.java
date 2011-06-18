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
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.RegionNavigatorDialog;
import org.broad.igv.ui.util.FileChooser;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.Collection;

/**
 * @author jrobinso
 */
public class RegionsBaseMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(RegionsBaseMenuAction.class);

    IGV mainFrame;

    protected enum Direction {

        IMPORT, EXPORT

    }

    ;

    public RegionsBaseMenuAction(String name, Icon icon, int mnemonic) {
        super(name, icon, mnemonic);
    }

    public final void importExportRegionsOfInterest(Direction direction) {
        //dhmay adding 20110505: There was an intermittent bug in which the regions list was not getting
        //synched when the user made changes to the regions table in the region navigator dialog.  I'm addressing
        //this with a change to RegionNavigatorDialog, but, due to the intermittent nature of the bug, I'm adding
        //a catchall here, too.  This synchs everything from the nav dialog.
        RegionNavigatorDialog navDialog = RegionNavigatorDialog.getActiveInstance();
        if (navDialog != null)
            navDialog.updateROIsFromRegionTable();

        File exportRegionDirectory = PreferenceManager.getInstance().getLastExportedRegionDirectory();
        if (exportRegionDirectory == null) {
            exportRegionDirectory = Globals.getUserDirectory();
        }

        FileChooser exportedRegionFileChooser = new FileChooser(exportRegionDirectory);
        String title = null;
        if (direction == Direction.EXPORT) {
            title = "Export Regions of Interest ...";
        } else {
            title = "Import Regions of Interest ...";
        }
        exportedRegionFileChooser.setDialogTitle(title);
        File file = selectExportedRegionsFile(exportedRegionFileChooser, new File("regions.bed"), direction == Direction.EXPORT);

        if (file == null) {
            return;
        }

        // Read or write the ROI file
        if (direction == Direction.EXPORT) {
            writeRegionsOfInterestFile(file);
        } else {
            readRegionsOfInterestFile(file);
        }
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
                            RegionOfInterest regionOfInterest = new RegionOfInterest(data[0],start, end, name);
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

    private File selectExportedRegionsFile(FileChooser exportedRegionFileChooser, File currentFile, boolean isSave) {

        exportedRegionFileChooser.setSelectedFile(currentFile);

        // Display the dialog
        if (isSave) {
            exportedRegionFileChooser.showSaveDialog(mainFrame.getMainFrame());
        } else {
            exportedRegionFileChooser.showOpenDialog(mainFrame.getMainFrame());
        }

        mainFrame.resetStatusMessage();
        File file = exportedRegionFileChooser.getSelectedFile();
        if (file != null) {
            File directory = exportedRegionFileChooser.getCurrentDirectory();
            if (directory != null) {
                PreferenceManager.getInstance().setLastExportedRegionDirectory(directory);
            }
        }

        return file;
    }

    private void writeRegionsOfInterestFile(File roiFile) {

        if (roiFile == null) {
            log.info("A blank Region of Interest export file was supplied!");
            return;
        }
        try {
            Collection<RegionOfInterest> regions = IGV.getInstance().getSession().getAllRegionsOfInterest();

            if (regions == null || regions.isEmpty()) {
                return;
            }

            // Create export file
            roiFile.createNewFile();
            PrintWriter writer = null;
            try {
                writer = new PrintWriter(roiFile);
                for (RegionOfInterest regionOfInterest : regions) {
                    Integer regionStart = regionOfInterest.getStart();
                    if (regionStart == null) {
                        // skip - null starts are bad regions of interest
                        continue;
                    }
                    Integer regionEnd = regionOfInterest.getEnd();
                    if (regionEnd == null) {
                        regionEnd = regionStart;
                    }

                    // Write info in BED format
                    writer.print(regionOfInterest.getChr());
                    writer.print("\t");
                    writer.print(regionStart);
                    writer.print("\t");
                    writer.print(regionEnd);

                    if (regionOfInterest.getDescription() != null) {
                        writer.print("\t");
                        writer.println(regionOfInterest.getDescription());
                    } else {
                        writer.println();
                    }
                }
            } finally {

                if (writer != null) {
                    writer.close();
                }
            }
        } catch (Exception e) {
            log.error("Failed to write Region of Interest export file!", e);
        }
    }
}
