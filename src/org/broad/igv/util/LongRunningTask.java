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
package org.broad.igv.util;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.WaitCursorManager.CursorToken;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.util.concurrent.*;

/**
 * Utility class for executing long running tasks in their own thread (i.e. not on the swing event thread).
 *
 * @author jrobinso
 */
public class LongRunningTask implements Callable {

    private static Logger log = Logger.getLogger(LongRunningTask.class);

    private static ExecutorService threadExecutor = Executors.newFixedThreadPool(10);
    private static ScheduledExecutorService schedule = Executors.newScheduledThreadPool(1);

    Runnable runnable;

    public static Future submit(Runnable runnable) {
        if (Globals.isBatch() || !SwingUtilities.isEventDispatchThread()) {
            runnable.run();
            return null;
        } else {

            return threadExecutor.submit(new LongRunningTask(runnable));
        }
    }

    public LongRunningTask(Runnable runnable) {
        this.runnable = runnable;
    }

    public Object call() throws Exception {

        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            runnable.run();
            return null;
        } catch (Exception e) {
            MessageUtils.showMessage("<html>Unexpected error: " + e.getMessage() + ".<br>See igv.log for more details");
            log.error("Exception running task", e);
            return null;
        } finally {
            //log.info("Removing wait cursor " + runnable.getName());
            WaitCursorManager.removeWaitCursor(token);
        }

    }
}

