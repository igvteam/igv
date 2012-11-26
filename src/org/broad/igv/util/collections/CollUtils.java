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
