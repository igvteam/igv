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
package org.broad.igv.ui;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.dev.db.DBManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.util.FileUtils;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collection;

/**
 * This thread is registered upon startup and will get executed upon exit.
 */
public class ShutdownThread extends Thread {

    private static Logger log = Logger.getLogger(ShutdownThread.class);


    public static void runS() {
        // Cleanup jnlp files
        if (Globals.IS_MAC) {
            FileUtils.cleanupJnlpFiles();
        }

        DBManager.shutdown();
        CommandListener.halt();
        cleanupBamIndexCache();
        if (IGV.hasInstance())
            IGV.getInstance().saveStateForExit();
    }

    @Override
    public void run() {
        runS();
    }

    private static void cleanupBamIndexCache() {
        File dir = DirectoryManager.getCacheDirectory();
        for (File f : dir.listFiles()) {
            f.delete();
        }
    }

    private static void writeRegionsOfInterestFile(File roiFile) {

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
