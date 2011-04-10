/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
package org.broad.igv.ui;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.RegionOfInterest;

import java.io.File;
import java.io.FileFilter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;

/**
 * This thread is registered upon startup and will get executed upon exit.
 */
public class ShutdownThread extends Thread {

    private static Logger log = Logger.getLogger(ShutdownThread.class);


    @Override
    public void run() {

        // Cleanup jnlp files
        if (Globals.IS_MAC) {
            File desktop = new File(System.getProperty("user.home") + "/Desktop");
            if (desktop.exists() && desktop.isDirectory()) {
                cleanupJnlpFiles(desktop);
            }
            File downloads = new File(System.getProperty("user.home") + "/Downloads");
            if (downloads.exists() && downloads.isDirectory()) {
                cleanupJnlpFiles(downloads);
            }

            if (IGV.getInstance().getSession().getAllRegionsOfInterest().size() > 0) {
                File file = new File("Regions_" + System.currentTimeMillis() + ".bed");

                if (file != null) {
                    writeRegionsOfInterestFile(file);
                }

            }


        }

    }


    /**
     * Cleanup extra jnlp files.  This method is written specifcally for Mac os.
     */
    public static void cleanupJnlpFiles(File desktop) {

        if (desktop.exists() && desktop.isDirectory()) {
            File[] jnlpFiles = desktop.listFiles(new FileFilter() {

                public boolean accept(File arg0) {
                    return arg0.getName().startsWith("igv") && arg0.getName().endsWith(".jnlp");
                }
            });

            // Sort files by ascending version number
            Arrays.sort(jnlpFiles, new Comparator<File>() {

                public int compare(File file1, File file2) {
                    return getVersionNumber(file1.getName()) - getVersionNumber(file2.getName());
                }

                private int getVersionNumber(String fn) {
                    int idx = fn.indexOf(".jnlp");
                    int idx2 = fn.lastIndexOf("-");
                    if (idx2 < 0) {
                        return 0;
                    } else {
                        try {
                            return Integer.parseInt(fn.substring(idx2 + 1, idx));

                        } catch (NumberFormatException numberFormatException) {
                            return 0;
                        }
                    }

                }
            });

            // Delete all but the highest version (newest) jnlp file
            for (int i = 0; i < jnlpFiles.length - 1; i++) {
                jnlpFiles[i].delete();
            }

            // Strip the version nuber fro the newest file
            if (jnlpFiles.length > 1) {
                File newestFile = jnlpFiles[jnlpFiles.length - 1];
                String fn = newestFile.getName();
                int dotIndex = fn.indexOf(".jnlp");
                int dashIndex = fn.lastIndexOf("-");
                if (dashIndex > 1) {
                    String newName = fn.substring(0, dashIndex) + fn.substring(dotIndex);
                    newestFile.renameTo(new File(newestFile.getParentFile(), newName));
                }
            }
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
