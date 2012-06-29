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
package org.broad.igv.util.collections;

import org.apache.log4j.Logger;

import java.lang.ref.SoftReference;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.WeakHashMap;

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

