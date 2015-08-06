/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Fred Hutchinson Cancer Research Center and Broad Institute
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
