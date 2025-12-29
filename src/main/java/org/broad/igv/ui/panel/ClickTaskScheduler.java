package org.broad.igv.ui.panel;

import org.broad.igv.ui.UIConstants;

import java.util.TimerTask;

/**
 *
 * A utility class for sceduling single-click actions "in the future",
 *
 * @author jrobinso
 * @date Dec 17, 2010
 */
public class ClickTaskScheduler {

    private TimerTask currentClickTask;

    public void cancelClickTask() {
        if (currentClickTask != null) {
            currentClickTask.cancel();
            currentClickTask = null;
        }
    }

    public void scheduleClickTask(TimerTask task) {
        cancelClickTask();
        currentClickTask = task;
        (new java.util.Timer()).schedule(currentClickTask, UIConstants.getDoubleClickInterval());
    }
}
