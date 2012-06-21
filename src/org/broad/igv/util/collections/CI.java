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

/**
 * Class for case-insensitive collections
 */
public class CI{

    private CI(){};

    private static String convertKey(String key){
        if(key == null){
            return null;
        }
        return key.toLowerCase();
    }

    public static class CIHashSet extends HashSet<String>{
        public CIHashSet(){
            super();
        }

        @Override
        public boolean add(String key){
            return super.add(convertKey(key));
        }

        public boolean remove(String key){
            return super.remove(convertKey(key));
        }

        public boolean contains(String key){
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

        public CIHashMap(){
            super();
        }


        @Override
        public V put(String key, V value){
            return super.put(convertKey(key), value);
        }

        public V remove(String key){
            return super.remove(convertKey(key));
        }

        public V get(String key){
            return super.get(convertKey(key));
        }

        public boolean containsKey(String key){
            return super.containsKey(convertKey(key));
        }

    }


}
