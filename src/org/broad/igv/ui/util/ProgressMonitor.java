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

package org.broad.igv.ui.util;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;


/**
 * TODO This class probably shouldn't exist, bet we could use built-in progress monitors instead
 * @author eflakes
 */
public class ProgressMonitor {

    final public static String PROGRESS_PROPERTY = "PROGRESS_PROPERTY";
    final public static String STATUS_PROPERTY = "STATUS_STATUS";

    private boolean isReady = false;
    private int oldValue = 0;
    private PropertyChangeSupport propertyChangeSupport;

    private String oldStatus = "";

    public ProgressMonitor() {
        propertyChangeSupport = new PropertyChangeSupport(this);
        setReady(true);
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        propertyChangeSupport.addPropertyChangeListener(listener);
    }

    /**
     * Sends an event to update the listening progress bar.
     *
     * @param value This value is the percent amount to add to the current
     *              progress in the progress bar (which goes from 0 to 100).
     */
    public void fireProgressChange(int value) {

        if (isReady) {
            final int newValue = oldValue + value;

            PropertyChangeEvent event =
                    new PropertyChangeEvent(
                            this,
                            PROGRESS_PROPERTY,
                            oldValue,
                            newValue);
            propertyChangeSupport.firePropertyChange(event);
            oldValue = newValue;
        }
    }

    public void updateStatus(String newStatus){
        if(isReady){
            PropertyChangeEvent event = new PropertyChangeEvent(this, STATUS_PROPERTY,
                    oldStatus, newStatus);
            propertyChangeSupport.firePropertyChange(event);
            oldStatus = newStatus;
        }
    }

    /**
     * Sets whether or not this class will respond to progress requests.
     *
     * @param ready
     */
    public void setReady(boolean ready) {
        isReady = ready;
    }

    /**
     * returns whether or not this class will respond to progress requests.
     *
     * @return
     */
    public boolean isReady() {
        return isReady;
    }
}