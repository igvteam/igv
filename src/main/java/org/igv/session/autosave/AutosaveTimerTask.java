package org.igv.session.autosave;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.session.Session;
import org.igv.ui.IGV;
import org.igv.ui.util.MessageUtils;

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
            // If autosaving fails, notify the user and stop autosaving
            final String message = "Failure while trying to autosave session. Timed autosave will be disabled until IGV is restarted.";
            log.error(message, err);
            MessageUtils.showMessage(message + "\n" + err);
            igv.stopTimedAutosave();
        }

    }
}
