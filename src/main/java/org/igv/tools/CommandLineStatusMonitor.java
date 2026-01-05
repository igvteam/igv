package org.igv.tools;

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
