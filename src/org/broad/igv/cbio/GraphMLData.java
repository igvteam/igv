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

/**
 * And individual datapoint in a GraphML file
 * We allow only those types allowed by GraphML:
 * boolean, int, long, double, float, String.
 * User: jacob
 * Date: 2012/02/01
 */
public class GraphMLData {

    //private String key;
    private String value;
    protected String type;
    //protected String defaultValue;

    private GraphMLData(Object value) {
        Class clazz = value.getClass();
        if (!clazz.isPrimitive() && !clazz.equals(String.class)) {
            throw new IllegalArgumentException("Value must be primitive or string");
        }
        //this.key = key;
        this.value = "" + value;
        this.type = clazz.getSimpleName().toLowerCase();
    }

    //We use these constructors for compile-time type safety
    //may get rid of them in the future
    public GraphMLData(String value) {
        this((Object) value);
    }

    public GraphMLData(long value) {
        this((Object) value);
    }

    public GraphMLData(int value) {
        this((Object) value);
    }

    public GraphMLData(double value) {
        this((Object) value);
    }

    public GraphMLData(float value) {
        this((Object) value);
    }

    public GraphMLData(boolean value) {
        this((Object) value);
    }

    public String getValue() {
        return this.value;
    }

    public String getType() {
        return this.type;
    }

}
