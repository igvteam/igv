package org.igv.ui.util;

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

    /**
     * Sends an event to update the listening progress bar.
     *
     * @param newValue This value is the percent amount completed (which goes from 0 to 100).
     */
    public void fireProgress(final int newValue) {

        if (isReady) {

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
