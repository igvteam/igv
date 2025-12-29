package org.igv.event;

import org.igv.logging.*;

import java.util.*;

/**
 * Ludicrously simple event bus -- its all we need.
 */
public class IGVEventBus {

    final Map<Class<?>, Set<IGVEventObserver>> observerMap;

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

    public synchronized void subscribe(Class<?> eventClass, IGVEventObserver observer) {
        Set<IGVEventObserver> observerSet = observerMap.computeIfAbsent(eventClass, k -> Collections.newSetFromMap(new WeakHashMap<>()));
        observerSet.add(observer);
    }

    /**
     * Unsubscribe observer from all observer lists.  Note that if this method is not called the observer should
     * still be eligible for garbage collection, but its good practice to call this nonetheless.
     */
    public synchronized void unsubscribe(IGVEventObserver observer) {
        for (Set<IGVEventObserver> observerSet : observerMap.values()) {
            observerSet.remove(observer);
        }
    }

    public void post(IGVEvent event) {
        Set<IGVEventObserver> observerSet = observerMap.get(event.getClass());
        if (observerSet != null) {
            // Make a copy in case original is modified during loop
            Collection<IGVEventObserver> observers = new ArrayList<>(observerSet);
            for (IGVEventObserver observer : observers) {
                observer.receiveEvent(event);
            }
        }
    }

}


