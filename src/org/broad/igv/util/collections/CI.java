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

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;

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
