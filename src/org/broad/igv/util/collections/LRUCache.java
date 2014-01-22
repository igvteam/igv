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

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * A least-recently-used cache with a maximum size that can be altered.
 * @author jrobinso
 */
public class LRUCache<K, V> {

    //Logger log = Logger.getLogger(LRUCache.class);

    private AtomicInteger maxEntries;

    private Map<K, V> map;


    public LRUCache(int max) {
        this.maxEntries = new AtomicInteger(max);
    }

    public void setMaxEntries(int max){
        this.maxEntries.set(max);
    }

    private void createMap() {
        map = Collections.synchronizedMap(
                new LinkedHashMap<K, V>(16, 0.75f, true) {
                    @Override
                    protected boolean removeEldestEntry(Map.Entry eldest) {
                        return (size() > maxEntries.get());
                    }
                });
    }

    private Map<K, V> getMap() {
        if (map == null) {
            createMap();
        }
        return map;
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

    public Set<K> keySet() {
        return getMap().keySet();
    }

    public Collection<V> values(){
        return getMap().values();
    }

    public void putAll(Map<K, V> map) {
        getMap().putAll(map);
    }

    public void putAll(LRUCache<K, V> intervals) {
        this.putAll(intervals.getMap());
    }
}

