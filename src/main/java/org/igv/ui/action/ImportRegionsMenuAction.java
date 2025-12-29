/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.action;

import org.igv.logging.*;
import org.igv.feature.RegionOfInterest;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.util.FileDialogUtils;
import org.igv.ui.util.MessageUtils;
import org.igv.ui.util.UIUtilities;

import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

/**
 * @author jrobinso
 */
public class ImportRegionsMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(ImportRegionsMenuAction.class);
    IGV igv;

    public ImportRegionsMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
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
            log.warn("No regions found in file");
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
                            String chr = data[0];
                            Genome genome = GenomeManager.getInstance().getCurrentGenome();
                            if(genome != null) {
                                chr = genome.getCanonicalChrName(chr);
                            }
                            RegionOfInterest regionOfInterest = new RegionOfInterest(chr, start, end, name);

                            igv.addRegionOfInterest(regionOfInterest);
                        } catch (NumberFormatException numberFormatException) {
                            log.error("Error importing regions of interest", numberFormatException);
                        }
                    }
                }
            } finally {

                if (reader != null) {
                    reader.close();
                }
                igv.repaint();
            }
        } catch (Exception e) {
            log.error("Failed to write Region of Interest export file!", e);
        }
    }

}
