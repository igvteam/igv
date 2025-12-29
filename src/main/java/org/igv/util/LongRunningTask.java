package org.igv.util;

import org.igv.logging.*;
import org.igv.Globals;
import org.igv.ui.IGV;
import org.igv.ui.WaitCursorManager;
import org.igv.ui.WaitCursorManager.CursorToken;
import org.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.util.concurrent.*;

/**
 * Utility class for executing long running tasks in their own thread (i.e. not on the swing event thread).
 *
 * @author jrobinso
 */
public class LongRunningTask implements Callable<Void> {

    private static final Logger log = LogManager.getLogger(LongRunningTask.class);

    private static final ExecutorService threadExecutor = Executors.newFixedThreadPool(5);

    Runnable runnable;

    public static Executor getThreadExecutor() {
        return threadExecutor;
    }

    public static Future<Void> submit(Runnable runnable) {
        if (Globals.isBatch()) {
            runnable.run();
            return null;
        } else {
            return threadExecutor.submit(new LongRunningTask(runnable));
        }
    }

    private LongRunningTask(Runnable runnable) {
        this.runnable = runnable;
    }

    @Override
    public Void call() {

        CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            runnable.run();
        } catch (Throwable e) {
            MessageUtils.showMessage("<html>Unexpected error: " + e.getMessage() + ".<br>See igv.log for more details");
            log.error("Exception running task", e);
        } finally {
            WaitCursorManager.removeWaitCursor(token);

            synchronized (IGV.getInstance()) {
                IGV.getInstance().notifyAll();
            }
        }
        return null;

    }


}

