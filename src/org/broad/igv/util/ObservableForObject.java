/**
 * Copyright (c) 2010-2011 by Fred Hutchinson Cancer Research Center.  All Rights Reserved.

 * This software is licensed under the terms of the GNU Lesser General
 * Public License (LGPL), Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.

 * THE SOFTWARE IS PROVIDED "AS IS." FRED HUTCHINSON CANCER RESEARCH CENTER MAKES NO
 * REPRESENTATIONS OR WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED,
 * INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS,
 * WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL FRED HUTCHINSON CANCER RESEARCH
 * CENTER OR ITS TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR
 * ANY DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR
 * CONSEQUENTIAL DAMAGES, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS,
 * REGARDLESS OF  WHETHER FRED HUTCHINSON CANCER RESEARCH CENTER SHALL BE ADVISED,
 * SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
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
