/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *  
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *  
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.cbio;

import java.util.*;

/**
 * Class for generating DataKeys used elements in a graph.
 * We store all keys generated for later retrieval when
 * creating schema.
 * <p/>
 * User: jacob
 * Date: 2012/02/02
 */
public class KeyFactory {


    private Map<String, DataKey> keyMap;
    private String _for;
    //Types get wrapped when cast to objects
    private static final Set<Class> allowedTypes = new HashSet<Class>(6);

    static {
        Class[] all_types = new Class[]{Long.class, Integer.class, Float.class, Double.class, Boolean.class, String.class};
        for (Class clz : all_types) {
            allowedTypes.add(clz);
        }
    }

    /**
     * @param _for "graph", "all", "node", "edge"
     */
    public KeyFactory(String _for) {
        this._for = _for;
        this.keyMap = new HashMap<String, DataKey>();
    }

    public String getFor() {
        return this._for;
    }

    public Collection<DataKey> getKeySet() {
        return this.keyMap.values();
    }

    public DataKey getDataKey(String name, Object value) {
        Class clazz = value.getClass();
        if (!allowedTypes.contains(clazz)) {
            throw new IllegalArgumentException("Value is type " + clazz +
                    " which is not allowed. Only long, int, float, double, boolean, and string are allowed.");
        }

        if (keyMap.containsKey(name)) {
            return keyMap.get(name);
        } else {
            DataKey key = new DataKey(name, clazz.getSimpleName());
            keyMap.put(name, key);
            return key;
        }
    }

}
