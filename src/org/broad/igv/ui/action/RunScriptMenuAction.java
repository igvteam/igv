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
import org.broad.igv.PreferenceManager;
import org.broad.igv.batch.BatchRunner;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.FileChooserDialog;
import org.broad.igv.ui.util.FileDialogUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;


public class RunScriptMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadFilesMenuAction.class);
    IGV mainFrame;

    public RunScriptMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.LOAD_TRACKS_TOOLTIP);
    }

    public void actionPerformed(ActionEvent e) {

        if (e.getActionCommand().equalsIgnoreCase("run batch script...")) {
            File script = chooseScriptFile();
            if (script != null && script.isFile()) {
                final BatchRunner bRun = new BatchRunner(script.getPath());

                SwingWorker worker = new SwingWorker() {

                    @Override
                    protected Object doInBackground() throws Exception {
                        bRun.run();
                        return null;
                    }
                };

                worker.execute();

            }
        }
    }


    private File chooseScriptFile() {

        File lastDirectoryFile = PreferenceManager.getInstance().getLastTrackDirectory();
        File scriptFile = FileDialogUtils.chooseFile("Select Script", lastDirectoryFile, FileDialog.LOAD);

        if (scriptFile != null) {
            // Store the last accessed file location
            PreferenceManager.getInstance().setLastTrackDirectory(scriptFile.getParentFile());
        }

        mainFrame.resetStatusMessage();
        return scriptFile;
    }
}
