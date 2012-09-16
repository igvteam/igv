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

import org.apache.commons.collections.Predicate;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * Utility methods which are useful for collections.
 * <p/>
 * User: jacob
 * Date: 2012-Aug-30
 */
public class CollUtils {


    /**
     * Filters the provided collection in place. Only those objects for which
     * predicate(object) returns true will be kept. Collection must support object
     * removal.
     *
     * @param objects
     * @param predicate
     * @param <T>
     */
    public static <T> void filter(Collection<? extends T> objects, Predicate<T> predicate) {
        if (objects == null) return;
        Iterator<? extends T> iter = objects.iterator();
        while (iter.hasNext()) {
            T next = iter.next();
            if (!predicate.evaluate(next)) {
                iter.remove();
            }
        }
    }

    /**
     * Filters the provided collection in place. Only those objects for which
     * predicate(object) returns true will be kept. Collection must support object
     * removal.
     *
     * @param objects
     * @param predicate
     * @param <T>
     */
    public static <T> Collection<T> filteredCopy(Collection<? extends T> objects, Predicate<T> predicate) {
        if (objects == null) return null;
        Collection<T> coll = new ArrayList<T>(objects.size());
        Iterator<? extends T> iter = objects.iterator();
        while (iter.hasNext()) {
            T next = iter.next();
            if (predicate.evaluate(next)) {
                coll.add(next);
            }
        }
        return coll;
    }


    /**
     * Move the element of the list from its current location
     * to the specified index. If {@code list} does not contain
     * {@code listItem}, has no effect.
     *
     * @param list
     * @param index
     * @param listItem
     * @param <T>
     * @return True for success
     */
    public static <T> boolean moveInList(List<T> list, int index, T listItem) {
        boolean moved = list.remove(listItem);
        if (moved) {
            list.add(index, listItem);
        }
        return moved;
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
