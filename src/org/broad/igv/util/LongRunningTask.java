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

import java.util.concurrent.*;

/**
 * @author jrobinso
 */
public class LongRunningTask implements Callable {

    private static Logger log = Logger.getLogger(LongRunningTask.class);

    private static ExecutorService threadExecutor = Executors.newFixedThreadPool(10);
    private static ScheduledExecutorService schedule = Executors.newScheduledThreadPool(1);

    NamedRunnable runnable;

    public static Future submit(NamedRunnable runnable) {
        if (Globals.isBatch() || !IGV.getInstance().isStartupComplete()) {
            log.debug("Running: " + runnable.getName());
            runnable.run();
            log.debug("Finished : " + runnable.getName());
            return null;
        } else {

            return threadExecutor.submit(new LongRunningTask(runnable));
        }
    }


    /**
     * Schedule a task for execution in the future.
     *
     * @param runnable
     * @param time
     * @return
     */
    public static Future schedule(NamedRunnable runnable, long time) {
        return schedule.schedule(new LongRunningTask(runnable), time, TimeUnit.MILLISECONDS);
    }

    public LongRunningTask(NamedRunnable runnable) {
        this.runnable = runnable;
    }

    public Object call() throws Exception {

        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            //log.info("Running " + runnable.getName());
            long t0 = System.currentTimeMillis();
            runnable.run();
            if (log.isDebugEnabled()) {
                long dt = System.currentTimeMillis() - t0;
                log.debug(runnable.getName() + "  time= " + dt);
            }
            return null;
        } catch(Exception e) {
            log.error("Exception running task: " + runnable.getName(), e);
            return null;
        }
        finally {
            //log.info("Removing wait cursor " + runnable.getName());
            WaitCursorManager.removeWaitCursor(token);
        }

    }
}

