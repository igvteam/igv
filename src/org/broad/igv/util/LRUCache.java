/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.util;

import org.apache.log4j.Logger;

import java.lang.ref.ReferenceQueue;
import java.lang.ref.SoftReference;
import java.lang.ref.WeakReference;
import java.util.*;

/**
 * An LRUCache using SoftReferences
 * <p/>
 * This implementation was based on ideas in the following article by Dr. Heinz Kabutz
 * <p/>
 * http://archive.devx.com/java/free/articles/Kabutz01/Kabutz01-1.asp
 *
 * @author jrobinso
 */
public class LRUCache<K, V> {

    Logger log = Logger.getLogger(LRUCache.class);

    private final int minEntries;

    private final int maxEntries;

    private LinkedHashMap<K, SoftValue<K, V>> theMap;

    /**
     * The FIFO list of hard references, order of last inserted.
     */
    private final LinkedList<V> hardCache = new LinkedList<V>();

    private final ReferenceQueue queue = new ReferenceQueue();

    public LRUCache(int max) {
        this(0, max);
    }

    public LRUCache(int min, int max) {
        this.minEntries = min;
        this.maxEntries = max;
        theMap = new LinkedHashMap<K, SoftValue<K, V>>() {
            @Override
            protected boolean removeEldestEntry(Map.Entry eldest) {
                return (size() > maxEntries);
            }
        };
    }

    public V put(K k, V v) {
        processQueue(); // throw out garbage collected values first
        SoftValue<K, V> softValue = theMap.put(k, new SoftValue(v, k, queue));

        if (minEntries > 0) {
            hardCache.addFirst(v);
            if (hardCache.size() > minEntries) {
                hardCache.removeLast();
            }
        }

        return softValue == null ? null : softValue.get();

    }

    public V get(Object key) {
        processQueue();
        SoftValue<K, V> softValue = theMap.get(key);
        if (softValue == null) {
            return null;
        } else {
            V value = softValue.get();
            if (value == null) {
                theMap.remove(key);
            }
            return value;
        }
    }

    public boolean containsKey(Object o) {
        return theMap.containsKey(o);
    }

    public boolean isEmpty() {
        return theMap.isEmpty();    //To change body of overridden methods use File | Settings | File Templates.
    }

    public void clear() {
        theMap.clear();
        while (queue.poll() != null) ;
    }

    /**
     * Here we go through the ReferenceQueue and remove garbage
     * collected SoftValue objects from the HashMap by looking them
     * up using the SoftValue.key data member.
     */
    private void processQueue() {

        SoftValue sv;
        while ((sv = (SoftValue) queue.poll()) != null) {
            theMap.remove(sv.key); // we can access private data!
            System.out.println("Removed " + sv.key);
        }

    }

    /**
     * We define our own subclass of SoftReference which contains
     * not only the value but also the key to make it easier to find
     * the entry in the HashMap after it's been garbage collected.
     */
    private static class SoftValue<K, V> extends SoftReference<V> {
        private final K key; // always make data member final

        /**
         * Did you know that an outer class can access private data
         * members and methods of an inner class?  I didn't know that!
         * I thought it was only the inner class who could access the
         * outer class's private information.  An outer class can also
         * access private members of an inner class inside its inner
         * class.
         */
        private SoftValue(V v, K key, ReferenceQueue q) {
            super(v, q);
            this.key = key;
        }
    }


}

