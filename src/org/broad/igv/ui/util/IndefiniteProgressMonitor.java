/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.util;

import java.util.Timer;
import java.util.TimerTask;

/**
 * @author jrobinso
 */
public class IndefiniteProgressMonitor extends ProgressMonitor {

    int cycleTime;
    Timer timer;

    private static final int DEFAULT_CYCLE_TIME = 60;

    public IndefiniteProgressMonitor(){
        this(DEFAULT_CYCLE_TIME);
    }

    private IndefiniteProgressMonitor(int cycleTime) {
        this.cycleTime = cycleTime;
        timer = new Timer();
        setReady(true);
    }

    public void start() {
        timer.schedule(new CycleTask(), 0, 1000);
    }

    public void stop() {
        timer.cancel();
        fireProgressChange(100);
    }

    class CycleTask extends TimerTask {

        boolean stop = false;
        int progress = 0;
        int progressIncrement = 0;
        int direction = 1;
        long lastTime = System.currentTimeMillis();

        @Override
        public void run() {


            UIUtilities.invokeOnEventThread(new Runnable() {
                public void run() {
                    fireProgressChange(progressIncrement);
                }
            });


            long t = System.currentTimeMillis();
            progressIncrement = (int) (direction * (t - lastTime) / (10 * cycleTime));
            progress += progressIncrement;
            if (progress >= 90) {
                progress = 99;
                direction = -1;
            } else if (progress <
                    0) {
                progress = 1;
                direction = 1;
            }
            lastTime = t;
        }
    }
}
