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