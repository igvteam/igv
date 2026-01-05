package org.igv.ui.action;

import org.igv.logging.*;
import org.igv.DirectoryManager;
import org.igv.feature.RegionOfInterest;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.IGV;
import org.igv.ui.panel.RegionNavigatorDialog;
import org.igv.ui.util.FileDialogUtils;
import org.igv.ui.util.UIUtilities;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.PrintWriter;
import java.util.Collection;

/**
 * @author jrobinso
 */
public class ExportRegionsMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(ExportRegionsMenuAction.class);
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

        File exportRegionDirectory = PreferencesManager.getPreferences().getLastExportedRegionDirectory();
        if (exportRegionDirectory == null) {
            exportRegionDirectory = DirectoryManager.getUserDefaultDirectory();
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
            //log.warn("A blank Region of Interest export file was supplied!");
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
