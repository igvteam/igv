package org.broad.igv.session;

import org.broad.igv.DirectoryManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;

import java.io.File;
import java.io.IOException;
import java.util.TimerTask;

/**
 * Triggers a session autosave when dictated by a timer
 */
public class AutosaveTimerTask extends TimerTask {

    private static Logger log = LogManager.getLogger(AutosaveTimerTask.class);

    private IGV igv;

    /**
     * Creates a new task which will save the current session for the specified IGV
     * instance when triggered by a timer
     * @param igv We'll autosave the session from this IGV instance
     */
    public AutosaveTimerTask(IGV igv) {
        this.igv = igv;
    }

    /**
     * Saves the current session for this listener's IGV instance
     */
    @Override
    public void run() {
        // Get the current session so we can save it
        Session session = igv.getSession();
        try {
            // Create a new autosave file
            File autosaveFile = DirectoryManager.getNewSessionAutosaveFile();
            // Save to the file
            (new SessionWriter()).saveSession(session, autosaveFile);
        }
        catch (IOException err) {
            log.error("Failed to autosave session", err);
        }

    }
}
