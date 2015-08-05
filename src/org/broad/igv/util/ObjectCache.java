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

/*
 * ImageCache.java
 *
 * Created on December 1, 2007, 9:59 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.util;

import org.broad.igv.sam.Alignment;

import java.lang.ref.SoftReference;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;

/**
 *
 */
public class ObjectCache<K, V> {

    /**
     * Maximum number of objects to cache
     */
    private int maxSize = 1000;
    private LinkedHashMap<K, SoftReference<V>> map = new LinkedHashMap<K, SoftReference<V>>();

    public ObjectCache() {
        this(50);
    }

    public ObjectCache(int maxSize) {
        this.maxSize = maxSize;
        map = new LinkedHashMap<K, SoftReference<V>>(maxSize);
    }

    public void put(K key, V image) {
        if (map.size() == maxSize) {
            // Map has reached maximum size.  Remove the first(oldest) entry.
            // 
            for (K k : map.keySet()) {
                map.remove(k);
                break;
            }
        }
        SoftReference<V> SoftReference = new SoftReference<V>(image);
        map.put(key, SoftReference);
    }

    public V get(K key) {

        V image = null;
        SoftReference<V> SoftReference = map.get(key);

        if (SoftReference != null) {
            image = SoftReference.get();
        }
        return image;
    }

    public Collection<K> getKeys() {
        return map.keySet();
    }

    public void remove(K key) {
        map.remove(key);
    }

    public boolean containsKey(K key) {
        return map.containsKey(key);
    }

    public void clear() {
        map.clear();
    }

    public int size() {
        return map.size();
    }

    public Iterator<K> iterator() {
        return map.keySet().iterator();
    }

    public int getMaxSize() {
        return maxSize;
    }

    public void setMaxSize(int maxSize) {
        this.maxSize = maxSize;
    }

}
