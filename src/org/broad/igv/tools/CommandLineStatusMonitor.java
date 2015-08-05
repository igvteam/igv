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
package org.broad.igv.tools;

import java.io.PrintStream;

/**
 * @author jrobinso
 */
public class CommandLineStatusMonitor implements StatusMonitor {

    private boolean interrupted = false;
    private double percentComplete = 0;

    private PrintStream outStream;

    public CommandLineStatusMonitor(PrintStream outStream){
        this.outStream = outStream;
    }

    public void setPercentComplete(double percentComplete) {
        this.percentComplete = Math.min(100, Math.max(0, percentComplete));

        if (percentComplete % 10 == 0) {
            outStream.println("" + percentComplete + "% ");
        } else {
            outStream.print(".");
        }
    }

    public void incrementPercentComplete(double increment) {
        this.percentComplete += increment;
        outStream.println("" + percentComplete + "%");
        if (percentComplete >= 100) {
            outStream.println("Done");
        }
    }

    /**
     * @return the interrupted
     */
    public boolean isInterrupted() {
        return interrupted;
    }

    /**
     * @param interrupted the interrupted to set
     */
    public void setInterrupted(boolean interrupted) {
        this.interrupted = interrupted;
    }
}
