
package org.igv.tools;

/**
 * @author jrobinso
 */
public interface StatusMonitor {

    public void setPercentComplete(double percentComplete);

    public void incrementPercentComplete(double increment);

    public boolean isInterrupted();

}
