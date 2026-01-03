package org.igv.util.collections;

import java.util.*;

/**
 * A Map implementation that treats String keys as case-insensitive.
 * Keys are converted to lowercase internally.
 *
 * @param <V> the type of mapped values
 */
public class CaseInsensitiveMap<V> implements Map<String, V> {

    private final Map<String, V> map = new HashMap<>();

    private static String lowerKey(Object key) {
        return key == null ? null : key.toString().toLowerCase();
    }

    @Override
    public int size() {
        return map.size();
    }

    @Override
    public boolean isEmpty() {
        return map.isEmpty();
    }

    @Override
    public boolean containsKey(Object key) {
        return map.containsKey(lowerKey(key));
    }

    @Override
    public boolean containsValue(Object value) {
        return map.containsValue(value);
    }

    @Override
    public V get(Object key) {
        return map.get(lowerKey(key));
    }

    @Override
    public V put(String key, V value) {
        return map.put(lowerKey(key), value);
    }

    @Override
    public V remove(Object key) {
        return map.remove(lowerKey(key));
    }

    @Override
    public void putAll(Map<? extends String, ? extends V> m) {
        for (Entry<? extends String, ? extends V> entry : m.entrySet()) {
            put(entry.getKey(), entry.getValue());
        }
    }

    @Override
    public void clear() {
        map.clear();
    }

    @Override
    public Set<String> keySet() {
        return map.keySet();
    }

    @Override
    public Collection<V> values() {
        return map.values();
    }

    @Override
    public Set<Entry<String, V>> entrySet() {
        return map.entrySet();
    }
}
