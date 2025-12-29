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
    private LinkedHashMap<K, SoftReference<V>> map;

    public ObjectCache() {
        this(50);
    }

    public ObjectCache(int maxSize) {
        this.maxSize = maxSize;
        map = new LinkedHashMap<>(maxSize);
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
