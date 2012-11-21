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

package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.batch.BatchRunner;
import org.broad.igv.ui.IGV;
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
