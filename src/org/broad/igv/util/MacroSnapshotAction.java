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
