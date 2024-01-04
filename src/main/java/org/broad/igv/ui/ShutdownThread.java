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
