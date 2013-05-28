/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
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
