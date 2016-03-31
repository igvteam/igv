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

    private static final ExecutorService threadExecutor = Executors.newFixedThreadPool(10);

    Runnable runnable;

    public static Executor getThreadExecutor() {
        return threadExecutor;
    }

    public static Future submit(Runnable runnable) {
        if (Globals.isBatch() || !SwingUtilities.isEventDispatchThread()) {
            runnable.run();
            return null;
        } else {
            return threadExecutor.submit(new LongRunningTask(runnable));
        }
    }

    private LongRunningTask(Runnable runnable) {
        this.runnable = runnable;
    }

    public Object call() throws Exception {

        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            runnable.run();
        } catch (Exception e) {
            MessageUtils.showMessage("<html>Unexpected error: " + e.getMessage() + ".<br>See igv.log for more details");
            log.error("Exception running task", e);
        } finally {
            //log.info("Removing wait cursor " + runnable.getName());
            WaitCursorManager.removeWaitCursor(token);

            synchronized (IGV.getInstance()) {
                IGV.getInstance().notifyAll();
            }
        }
        return null;

    }


}

