package org.broad.igv.session.autosave;

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.session.Session;
import org.broad.igv.ui.IGV;

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
            // Save the session to a new file in the autosave directory
            SessionAutosaveManager.saveTimedSessionAutosaveFile(session);
        }
        catch (IOException err) {
            log.error("Failed to autosave session", err);
        }

    }
}
