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
