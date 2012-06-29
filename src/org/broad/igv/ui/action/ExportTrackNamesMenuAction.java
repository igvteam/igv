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

package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

/**
 * @author jrobinso
 * @date Nov 7, 2010
 */
public class ExportTrackNamesMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(ClearRegionsMenuAction.class);
    IGV igv;

    public ExportTrackNamesMenuAction(String label, IGV mainFrame) {
        super(label, null);
        this.igv = mainFrame;
        setToolTipText(UIConstants.EXPORT_REGION_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        FileDialog fd = new FileDialog(igv.getMainFrame());
        fd.setModal(true);
        fd.setMode(FileDialog.SAVE);
        fd.setVisible(true);

        String fname = fd.getFile();
        if (fname == null) {
            return;
        }

        File outputFile = new File(fd.getDirectory(), fname);

        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new FileWriter(outputFile));

            final List<ReferenceFrame> referenceFrames = FrameManager.getFrames();
            if (referenceFrames.size() > 1) {
                pw.print("Sample");
                for (ReferenceFrame frame : referenceFrames) {
                    pw.print("\t" + frame.getName());
                }
                pw.println();
            }

            for (Track t : igv.getAllTracks()) {
                if (t.getTrackType() == TrackType.COPY_NUMBER || t.getTrackType() == TrackType.CNV) {
                    pw.print(t.getName());
                    for (ReferenceFrame frame : referenceFrames) {
                        //track.getRegionScore(chr, start, end, zoom, type, frame));
                        String chr = frame.getChrName();
                        int start = (int) frame.getOrigin();
                        int end = (int) frame.getEnd();
                        float score = t.getRegionScore(chr, start, end, frame.getZoom(), RegionScoreType.SCORE, frame.getName());
                        pw.print("\t" + score);
                    }
                    pw.println();
                }
            }
        } catch (IOException ex) {
            MessageUtils.showMessage("IO Error: " + ex.getMessage());
        } finally {
            if (pw != null) {
                pw.close();
            }
        }

    }
}
