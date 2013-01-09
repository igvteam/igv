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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.util;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;

import java.io.File;
import java.util.List;

/**
 * @deprecated What is the purpose of this class? -Jacob S
 * @author jrobinso
 */
@Deprecated
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
//        IGV mainFrame = IGV.getInstance();
//
//        if (filename == null) {
//            String locus = FrameManager.getDefaultFrame().getFormattedLocusString();
//            filename = locus.replaceAll(":", "_").replace("-", "_") + ".png";
//        }
//
//        // Repaint
//        IGV.getInstance().repaintDataAndHeaderPanels();
//        IGV.getInstance().repaintStatusAndZoomSlider();
//
//
//        File file = new File(outputDirectory, filename);
//        log.info("Snapshot: " + filename);
//        try {
//            mainFrame.createSnapshotNonInteractive(file);
//        } catch (IOException e) {
//            log.error(e);
//            throw new RuntimeException(e);
//        }
    }
}
