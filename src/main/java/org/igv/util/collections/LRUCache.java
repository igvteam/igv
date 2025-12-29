package org.igv.util.collections;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * A least-recently-used cache with a maximum size that can be altered.
 *
 * @author jrobinso
 */
public class LRUCache<K, V> {

    //Logger log = LogManager.getLogger(LRUCache.class);

    private AtomicInteger maxEntries;

    private Map<K, V> map;


    public LRUCache(int max) {
        this.maxEntries = new AtomicInteger(max);
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

    public Collection<V> values() {
        return getMap().values();
    }

    public void putAll(Map<K, V> map) {
        getMap().putAll(map);
    }

    public void putAll(LRUCache<K, V> intervals) {
        this.putAll(intervals.getMap());
    }

    /**
     * Test for memory leaks
     */

    public static void main(String[] args) {

        LRUCache<String, int [] > cache = new LRUCache<String, int []>(10);


        for (int i = 0; i < Integer.MAX_VALUE; i++) {

            cache.put(String.valueOf(i), new int[1000]);

            if(i % 1000000 == 0) {
                // Get the Java runtime
                Runtime runtime = Runtime.getRuntime();
                // Run the garbage collector
                runtime.gc();
                // Calculate the used memory
                long memory = runtime.totalMemory() - runtime.freeMemory();
                System.out.println("Used memory (kb): " + memory / 1000);
            }

        }
    }
}

