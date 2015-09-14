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

import com.google.common.base.Predicate;
import com.google.common.collect.Collections2;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Utility methods which are useful for collections.
 * <p/>
 * User: jacob
 * Date: 2012-Aug-30
 */
public class CollUtils {

    public static <T> List<T> filter(Collection<T> unfiltered, Predicate<? super T> predicate) {
        if (unfiltered == null) return null;
        Collection<T> filteredColl = Collections2.filter(unfiltered, predicate);
        List<T> filteredList = new ArrayList<T>(filteredColl);
        return filteredList;
    }

    /**
     * @param enumType
     * @param name
     * @return The Enum type with this exact {@code name}, or {@code defaultValue} if not found
     * @see #findValueOf(Class, String)
     */
    public static <T extends Enum<T>> T valueOf(Class<T> enumType, String name, T defaultValue) {
        try {
            return Enum.<T>valueOf(enumType, name);
        } catch (Exception e) {
            return defaultValue;
        }
    }

    /**
     * Search through an to find an object of matching value.
     * Will match either the text name of the enum, or the human readable name
     * Some our enums have a String value, intended to be human readable
     * <p/>
     * Comparisons with the enum name are case sensitive, comparisons
     * with the value (human readable name) are case insensitive
     *
     * @param collectionClass
     * @param value
     * @param <T>
     * @return
     */
    public static <T extends Enum<T> & Valued> T findValueOf(Class<T> collectionClass, String value) {
        if (value == null) {
            return null;
        }

        T[] enumConstants = collectionClass.getEnumConstants();
        if (enumConstants == null) {
            throw new IllegalArgumentException("Input must be an enum: " + collectionClass);
        }

        for (T valued : enumConstants) {
            if (value.equals(valued.name()) || value.equalsIgnoreCase(valued.getValue())) {
                return valued;
            }
        }

        return null;
    }

    public static interface Valued {
        String getValue();
    }
}
