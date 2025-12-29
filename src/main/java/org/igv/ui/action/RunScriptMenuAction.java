package org.igv.ui.action;

import org.igv.logging.*;
import org.igv.batch.BatchRunner;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.IGV;
import org.igv.ui.util.FileDialogUtils;
import org.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;


public class RunScriptMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(LoadFilesMenuAction.class);
    IGV igv;

    public RunScriptMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    /**
     * Run the batch script.  This is PURPOSELY run on the event dispatch thread to maintain absolute synchronization
     * @param e
     */
    public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equalsIgnoreCase("run batch script...")) {
            File [] script = chooseScriptFile();
            if (script != null && script.length > 0) {
                for(File f : script) {
                    final BatchRunner bRun = new BatchRunner(f.getPath(), igv);
                    try {
                        bRun.run();
                    } catch (Exception ex) {
                        MessageUtils.showMessage(Level.ERROR, "Error running batch script: " + ex.getMessage());
                    }
                }
            }
        }
    }


    private File [] chooseScriptFile() {

        File lastDirectoryFile = PreferencesManager.getPreferences().getLastTrackDirectory();
        File [] scriptFile = FileDialogUtils.chooseMultiple("Select Script", lastDirectoryFile, null);

        if (scriptFile != null && scriptFile.length > 0) {
            // Store the last accessed file location
            PreferencesManager.getPreferences().setLastTrackDirectory(scriptFile[0].getParentFile());
        }

        igv.resetStatusMessage();
        return scriptFile;
    }
}
