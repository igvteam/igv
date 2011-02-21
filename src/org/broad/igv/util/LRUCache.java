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
import java.util.*;

/**
 * @author jrobinso
 */
public class LRUCache<K, V> {

    Logger log = Logger.getLogger(LRUCache.class);



    private static Map<Object, LRUCache> instances = Collections.synchronizedMap(new WeakHashMap<Object, LRUCache>());

    public static void clearCaches() {
        for (LRUCache cache : instances.values()) {
            if (cache != null) {
                cache.clear();
            }
        }
    }



    private final int maxEntries;

    private SoftReference<Map<K, V>> mapReference;


    public LRUCache(Object source, int max) {
        instances.put(source, this);
        this.maxEntries = max;
        createMap();
    }

    private void createMap() {
        mapReference = new SoftReference(Collections.synchronizedMap(
                new LinkedHashMap<K, V>() {
                    @Override
                    protected boolean removeEldestEntry(Map.Entry eldest) {
                        return (size() > maxEntries);
                    }
                }));
    }

    private Map<K, V> getMap() {
        if (mapReference == null || mapReference.get() == null) {
            System.out.println("Creating map");
            createMap();
        }
        return mapReference.get();
    }

    public V put(K k, V v) {
        return getMap().put(k, v);
    }

    public V get(Object key) {
        return getMap().get(key);
    }

    public boolean containsKey(Object o) {
        return getMap().containsKey(o);
    }

    public boolean isEmpty() {
        return getMap().isEmpty();    //To change body of overridden methods use File | Settings | File Templates.
    }

    public void clear() {
        getMap().clear();
    }




}

