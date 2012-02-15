/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
