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
package org.broad.igv.util;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * @author jrobinso
 */
public class MacroSnapshotAction {

    public static File OUTPUT_DIRECTORY = new File(".");
    private static Logger log = Logger.getLogger(MacroSnapshotAction.class);

    /**
     * Loop through a list of loci creating a screenshot for each.  Method
     * assumes that data has been loaded.
     */
    final public static void doScreenshots(File regionFile, File outputDirectory) {
        doSnapshots(regionFile, outputDirectory);
    }

    final public static void doSnapshots(File regionFile, File outputDirectory) {

        final List<String> loci = ParsingUtils.loadRegions(regionFile);

        for (String locus : loci) {
            IGV.getInstance().goToLocus(locus);
            createSnapshot(outputDirectory, null);
        }
    }

    final public static synchronized void setOutputDirectory(String dir) {
        OUTPUT_DIRECTORY = new File(dir);
        if (!OUTPUT_DIRECTORY.exists()) {
            log.error("Warning: non existent directory: " + dir);
        }
    }


    final public static void doSnapshot(String filename) {
        createSnapshot(OUTPUT_DIRECTORY, filename);
    }

    private static synchronized void createSnapshot(File outputDirectory, String filename) {
        IGV mainFrame = IGV.getInstance();

        if (filename == null) {
            String locus = FrameManager.getDefaultFrame().getFormattedLocusString();
            filename = locus.replaceAll(":", "_").replace("-", "_") + ".png";
        }

        // Repaint
        IGV.getInstance().repaintDataAndHeaderPanels();
        IGV.getInstance().repaintStatusAndZoomSlider();


        File file = new File(outputDirectory, filename);
        log.info("Snapshot: " + filename);
        try {
            mainFrame.createSnapshotNonInteractive(file);
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeException(e);
        }
    }
}
