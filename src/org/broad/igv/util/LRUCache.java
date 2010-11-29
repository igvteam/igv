/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
package org.broad.igv.util;

import org.apache.log4j.Logger;

import java.lang.ref.WeakReference;
import java.util.*;

/**
 * @author jrobinso
 */
public class LRUCache<K, V> extends LinkedHashMap<K, V> {

    Logger log = Logger.getLogger(LRUCache.class);

    private int maxEntries = 100;

    private static Set<WeakReference<LRUCache>> instances = new HashSet();


    public LRUCache(int maxEntries) {
        this.maxEntries = maxEntries;
        instances.add(new WeakReference(this));
    }

    public static void clearCaches() {
        for (WeakReference<LRUCache> cache : instances) {
            LRUCache instance = cache.get();
            if (instance != null) {
                instance.clear();
            }
        }
        instances.clear();
    }

    public static void removeAllOldestEntries() {
        for (WeakReference<LRUCache> cache : instances) {
            LRUCache instance = cache.get();
            if (instance != null) {
                instance.removeOldestEntries();
            }
        }
    }

    @Override
    protected boolean removeEldestEntry(Map.Entry eldest) {
        if (size() > maxEntries) {
            if (log.isDebugEnabled()) {
                log.debug("Exceed cache capacity.");
            }
            return true;
        } else if (RuntimeUtils.getAvailableMemoryFraction() < 0.3) {
            log.info("Memory low.  Free cache entry");
            removeAllOldestEntries();

            // Per the LinkedHashMap contract return false, to prevent further modifications
            return false;
        } else {
            return false;
        }
    }

    private void removeOldestEntries() {
        // Remove oldest entries,  keeping the lower of the newest 5 or 1/2 of the entries.
        List<K> keys = new ArrayList(keySet());
        int midPoint = Math.min(5, keys.size() / 2);
        for (int i = midPoint; i < keys.size(); i++) {
            remove(keys.get(i));
        }
        System.gc();
    }
}

