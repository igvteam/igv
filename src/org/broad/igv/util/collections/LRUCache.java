/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
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

