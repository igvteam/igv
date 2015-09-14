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

/**
 * A "map-like" class that supports multiple values for a given key.
 * Note this violates the map contract,  thus "map-like".
 * Created to support GFF column 9 and Genbank tag-value pairs.
 */
public class MultiMap<K, V> {

    LinkedHashMap<K, Object> map;
    int size = 0;

    public MultiMap() {
        this(10);
    }

    public MultiMap(int size) {
        map = new LinkedHashMap<K, Object>(size);
    }

    public void put(K key, V value) {
        Object currentValue = map.get(key);
        if (currentValue == null) {
            map.put(key, value);
        } else if (currentValue instanceof List) {
            ((List) currentValue).add(value);
        } else {
            List<V> valueList = new ArrayList<V>();
            valueList.add((V) currentValue);
            valueList.add(value);
            map.put(key, valueList);
        }
        size++;
    }

    public V get(K key) {
        Object value = map.get(key);
        if (value == null) {
            return null;
        } else if (value instanceof List) {
            return ((List<V>) value).get(0);
        } else {
            return (V) value;
        }
    }

    public V remove(K key) {
        Object value = map.get(key);
        if (value == null) return null;
        if (value instanceof List) {
            size -= ((List) value).size();
            map.remove(key);
            return (V) ((List) value).get(0);
        } else {
            size--;
            map.remove(key);
            return (V) value;
        }
    }

    public int size() {
        return size;
    }

    public boolean containsKey(K key) {
        return map.containsKey(key);
    }

    public Set<K> keys() {
        return map.keySet();
    }

    private static final int MAX_CHARS_PER_LINE = 50;

    public void printHtml(StringBuffer buffer, int max) {

        if (map == null || map.isEmpty()) return;

        int count = 0;
        buffer.append("<br>");
        for (Map.Entry<K, Object> entry : map.entrySet()) {
            Object value = entry.getValue();
            if (value instanceof List) {
                for (V v : ((List<V>) value)) {
                    buffer.append("<b>" + entry.getKey() + "</b>");
                    buffer.append(":&nbsp;");
                    buffer.append(v.toString());
                    buffer.append("<br/>");
                    count++;
                }
            } else {
                buffer.append("<b>" + entry.getKey() + "</b>");
                buffer.append(":&nbsp;");
                String ts = lineWrapString(value.toString(), MAX_CHARS_PER_LINE);

                buffer.append(ts);
                buffer.append("<br/>");
                count++;
            }
            if (++count > max) {
                buffer.append("...");
                break;
            }

        }
    }

    private String lineWrapString(String input, int maxCharsPerLine) {
        int lines = input.length() / maxCharsPerLine + 1;
        if (lines == 1) return input;

        String result = input.substring(0, maxCharsPerLine);
        for (int lineNum = 1; lineNum < lines; lineNum++) {
            int start = lineNum * maxCharsPerLine;
            int end = Math.min(start + maxCharsPerLine, input.length());
            result += "<br/>" + input.substring(start, end);
        }
        return result;
    }

    public List<V> values() {
        List<V> allValues = new ArrayList<V>(map.size());
        for (Map.Entry<K, Object> entry : map.entrySet()) {
            Object value = entry.getValue();
            if (value instanceof List) {
                allValues.addAll((List) value);

            } else {
                allValues.add((V) value);
            }

        }
        return allValues;
    }

    public void addAll(MultiMap<K, V> attributes) {
        for (K key : attributes.keys()) {
            map.put(key, attributes.get(key));
        }
    }
}
