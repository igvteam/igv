/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
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

package org.broad.igv.event;


import org.apache.log4j.Logger;

import java.util.*;

/**
 * Ludicrously simple event bus -- its all we need.
 */
public class IGVEventBus {

    static final Logger log = Logger.getLogger(IGVEventBus.class);

    Map<Class, Set<IGVEventObserver>> observerMap;

    private static IGVEventBus instance;

    // This is not a singleton,  "instance" is the default bus.  Sashimi plot has its own bus.
    public static synchronized IGVEventBus getInstance() {
        if (instance == null) {
            instance = new IGVEventBus();
        }
        return instance;
    }

    public IGVEventBus() {
        this.observerMap = new HashMap<>();
    }

    public synchronized void subscribe(Class eventClass, IGVEventObserver observer) {

        Set<IGVEventObserver> observerSet = observerMap.get(eventClass);
        if (observerSet == null) {
            observerSet = Collections.newSetFromMap(new WeakHashMap<IGVEventObserver, Boolean>());
            observerMap.put(eventClass, observerSet);
        }
        observerSet.add(observer);
    }

    /**
     * Unsubscribe observer from all observer lists.  Not that if this method is not called the observer should
     * still be eligible for garbage collection, but its good practice to call this nonetheless.
     */
    public synchronized void unsubscribe(IGVEventObserver observer) {
        for (Set<IGVEventObserver> observerSet : observerMap.values()) {
            observerSet.remove(observer);
        }
    }

    public void post(Object event) {
        Set<IGVEventObserver> observerSet = observerMap.get(event.getClass());  // Make a copy in case original is modified during loop
        if (observerSet == null) {
            log.info("No observers for event type: " + event.getClass());
        } else {
            // Make a copy in case original is modified during loop
            Collection<IGVEventObserver> observers = new ArrayList<>(observerSet);
            for (IGVEventObserver observer : observers) {
                observer.receiveEvent(event);
            }
        }
    }

}


