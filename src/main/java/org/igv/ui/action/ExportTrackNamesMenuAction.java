package org.igv.ui.action;

import org.igv.logging.*;
import org.igv.track.RegionScoreType;
import org.igv.track.Track;
import org.igv.track.TrackType;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.MessageUtils;

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

    static Logger log = LogManager.getLogger(ClearRegionsMenuAction.class);
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
