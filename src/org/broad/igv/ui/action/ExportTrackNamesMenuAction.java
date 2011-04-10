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

package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.MessageUtils;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author jrobinso
 * @date Nov 7, 2010
 */
public class ExportTrackNamesMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(ClearRegionsMenuAction.class);
    IGV mainFrame;

    public ExportTrackNamesMenuAction(String label, IGV mainFrame) {
        super(label, null);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.EXPORT_REGION_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        FileDialog fd = new FileDialog(mainFrame.getMainFrame());
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

            for (Track t : mainFrame.getTrackManager().getAllTracks(false)) {
                pw.println(t.getName());
            }
        }
        catch (IOException ex) {
            MessageUtils.showMessage("IO Error: " + ex.getMessage());
        }
        finally {
            if (pw != null) {
                pw.close();
            }
        }

    }
}
