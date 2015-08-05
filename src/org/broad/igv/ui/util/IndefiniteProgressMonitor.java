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
