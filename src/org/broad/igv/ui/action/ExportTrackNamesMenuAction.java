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
