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
import org.broad.igv.DirectoryManager;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.RegionNavigatorDialog;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.UIUtilities;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.PrintWriter;
import java.util.Collection;

/**
 * @author jrobinso
 */
public class ExportRegionsMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(ExportRegionsMenuAction.class);
    IGV mainFrame;

    public ExportRegionsMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                exportRegionsOfInterest();
            }
        });
    }

    public final void exportRegionsOfInterest() {
        //dhmay adding 20110505: There was an intermittent bug in which the regions list was not getting
        //synched when the user made changes to the regions table in the region navigator dialog.  I'm addressing
        //this with a change to RegionNavigatorDialog, but, due to the intermittent nature of the bug, I'm adding
        //a catchall here, too.  This synchs everything from the nav dialog.
        RegionNavigatorDialog navDialog = RegionNavigatorDialog.getInstance();
        if (navDialog != null) navDialog.updateROIsFromRegionTable();

        File exportRegionDirectory = PreferenceManager.getInstance().getLastExportedRegionDirectory();
        if (exportRegionDirectory == null) {
            exportRegionDirectory = DirectoryManager.getUserDirectory();
        }

        String title = "Export Regions of Interest ...";
        File file = FileDialogUtils.chooseFile(title, exportRegionDirectory, new File("regions.bed"), FileDialog.SAVE);

        if (file == null) {
            return;
        }
        writeRegionsOfInterestFile(file);

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
