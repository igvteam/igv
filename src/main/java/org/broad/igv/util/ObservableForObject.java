package org.broad.igv.util;

import java.util.Observable;

/**
 * Contains an object of some type, and provides an accessor and a way of notifying observers that the thing
 * has changed.
 * <p/>
 * The only reason this class is necessary is to expose a way to get at setChanged(),
 * which is protected in Observable.  The typing is a convenience, and so is the storage of
 * the object that we really care about.
 */
public class ObservableForObject<T> extends Observable {

    protected T observedThing;

    public ObservableForObject(T observedList) {
        this.observedThing = observedList;
    }

    /**
     * Get the thing that is observed
     *
     * @return
     */
    public T getThing() {
        return observedThing;
    }

    /**
     * Change the stored object that is observed
     *
     * @param observedThing
     */
    public void setThing(T observedThing) {
        this.observedThing = observedThing;
    }

    public void clearChanged() {
        super.clearChanged();
    }

    public void setChanged() {
        super.setChanged();
    }

    /**
     * indicate that the object is changed and notify observers of the change
     */
    public void setChangedAndNotify() {
        setChanged();
        notifyObservers();
    }
}
