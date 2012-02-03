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

import java.util.HashMap;

/**
 * Named element to be used in a graph. If not null,
 * we use the name field as a unique identifier.
 * User: jacob
 * Date: 2012/02/01
 */

public abstract class BaseElement extends HashMap<DataKey, String> {

    private KeyFactory factory;
    private String name;

    public BaseElement(String name, KeyFactory factory) {
        super();
        this.name = name;
        this.factory = factory;
    }

    /**
     * Add data to this map. Note that we generate
     * a DataKey using the KeyFactory.
     *
     * @param name
     * @param value
     * @return
     */
    public String put(String name, Object value) {
        return this.put(this.factory.getDataKey(name, value), "" + value);
    }

    /**
     * Check if Edge contains a key with key named {@code name}
     *
     * @param name
     * @return
     */
    public boolean containsKey(String name) {
        //Don't use factory because we don't want to store the key
        DataKey key = new DataKey(name, "");
        return this.containsKey(key);
    }

    /**
     * Get the value with key named {@code name}, null if not found.
     *
     * @param name
     * @return
     */
    public String get(String name) {
        DataKey key = new DataKey(name, "");
        return this.get(key);
    }

    @Override
    public String remove(Object object) {
        throw new UnsupportedOperationException("Cannot remove data from edge");
    }

    @Override
    public boolean equals(Object o1) {
        if (this.name == null) {
            return super.equals(o1);
        }
        try {
            BaseElement okey = (BaseElement) o1;
            return this.name.equals(okey.getName());
        } catch (ClassCastException e) {
            return false;
        }
    }

    @Override
    public int hashCode() {
        if (this.name == null) {
            return super.hashCode();
        }
        return this.name.hashCode();
    }


    public String getName() {
        return name;
    }
}