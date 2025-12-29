/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.tools;

/**
 * @author jrobinso
 */
public interface StatusMonitor {

    public void setPercentComplete(double percentComplete);

    public void incrementPercentComplete(double increment);

    public boolean isInterrupted();

}
