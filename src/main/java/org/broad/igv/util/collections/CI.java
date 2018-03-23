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

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Class for case-insensitive collections. Needless to say,
 * the keys must be Strings.
 */
public class CI {

    private CI() {
    }

    ;

    private static String convertKey(Object key) {
        if (key == null) {
            return null;
        }
        return ((String) key).toLowerCase();
    }

    public static class CIHashSet extends HashSet<String> {
        public CIHashSet() {
            super();
        }

        @Override
        public boolean add(String key) {
            return super.add(convertKey(key));
        }

        @Override
        public boolean remove(Object key) {
            return super.remove(convertKey(key));
        }

        @Override
        public boolean contains(Object key) {
            return super.contains(convertKey(key));
        }
    }

    /**
     * A hashmap which uses strings for keys, that are
     * case insensitive. In practice this means all keys are
     * set to lowercase; the iterator will reflect that.
     * User: Jacob
     * Date: 2012-Jun-21
     */
    public static class CIHashMap<V> extends HashMap<String, V> {

        public CIHashMap() {
            super();
        }

        public CIHashMap(Map<String, V> inMap){
            for(Map.Entry<String, V> entry: inMap.entrySet()){
                this.put(entry.getKey(), entry.getValue());
            }
        }

        @Override
        public V put(String key, V value) {
            return super.put(convertKey(key), value);
        }

        @Override
        public V remove(Object key) {
            return super.remove(convertKey(key));
        }

        @Override
        public V get(Object key) {
            return super.get(convertKey(key));
        }

        @Override
        public boolean containsKey(Object key) {
            return super.containsKey(convertKey(key));
        }

    }

    public static class CILinkedHashMap<V> extends LinkedHashMap<String, V> {

        public CILinkedHashMap() {
            super();
        }

        public CILinkedHashMap(int size) {
            super(size);
        }


        @Override
        public V put(String key, V value) {
            return super.put(convertKey(key), value);
        }

        @Override
        public V remove(Object key) {
            return super.remove(convertKey(key));
        }

        @Override
        public V get(Object key) {
            return super.get(convertKey(key));
        }

        @Override
        public boolean containsKey(Object key) {
            return super.containsKey(convertKey(key));
        }

    }


}
