package org.igv.ui.action;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.IGV;
import org.igv.ui.util.FileDialogUtils;
import org.igv.util.LongRunningTask;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;

/**
 * This menu action classes is used for both "Open Session ..." and load recent
 * session menu items.  In the "Open Session..." the user has to specify a
 * session file through the file menu.  For "load recent"  the action is
 * instantiated with a specific session file.
 *
 * @author jrobinso
 */
public class OpenSessionMenuAction extends MenuAction {

    private static Logger log = LogManager.getLogger(OpenSessionMenuAction.class);
    private IGV igv;
    private String sessionFile = null;
    private boolean autoload = false;

    /**
     * Constructor for hardcoded session list at bottom of menu.   This convention is odd, legacy thing.
     * @param sessionFile
     * @param igv
     */
    public OpenSessionMenuAction(String sessionFile, IGV igv) {
        super(sessionFile);
        this.sessionFile = sessionFile;
        this.igv = igv;
        this.autoload = true;
    }

    public OpenSessionMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        if (sessionFile == null || autoload == false) {
            sessionFile = pickSessionFile();
        }
        if (sessionFile != null) {
            LongRunningTask.submit(() -> this.igv.loadSession(sessionFile, null));
        }
    }

    private static String pickSessionFile() {
        File lastSessionDirectory = PreferencesManager.getPreferences().getLastTrackDirectory();
        File tmpFile = FileDialogUtils.chooseFile("Open Session", lastSessionDirectory, JFileChooser.FILES_ONLY);

        final String result;
        if (tmpFile == null) {
            result = null;
        } else {
            result = tmpFile.getAbsolutePath();
            PreferencesManager.getPreferences().setLastTrackDirectory(tmpFile.getParentFile());
        }
        return result;
    }
}

