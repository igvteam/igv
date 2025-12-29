package org.broad.igv.ui;

import org.broad.igv.logging.*;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.Track;

/**
 * This thread is registered upon startup and will get executed upon exit.
 */
public class ShutdownThread extends Thread {

    private static Logger log = LogManager.getLogger(ShutdownThread.class);
    private static long oneDayMS = 24 * 60 * 60 * 1000;

    public static void runS() {
        log.info("Shutting down");
        CommandListener.halt();
        if (IGV.hasInstance()) {
            try {
                IGV.getInstance().saveStateForExit();
            } catch (Exception e) {
                log.error("Error saving session ", e);
            }
            PreferencesManager.getPreferences().setApplicationFrameBounds(IGV.getInstance().getMainFrame().getBounds());
            for (Track t : IGV.getInstance().getAllTracks()) {
                t.unload();
            }
        }
        LogFileHandler.getInstance().closeHandler();
    }

    @Override
    public void run() {
        runS();
    }

}
