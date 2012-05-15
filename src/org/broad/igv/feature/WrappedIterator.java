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

package org.broad.igv.feature;

import org.broad.tribble.CloseableTribbleIterator;

import java.util.Iterator;

/**
 * User: jacob
 * Date: 2012/05/15
 */
public class WrappedIterator<Feature> implements CloseableTribbleIterator {

    private Iterator<Feature> iterator;

    public WrappedIterator(Iterator<Feature> iterator) {
        this.iterator = iterator;
    }

    public void close() {

    }

    public Iterator<Feature> iterator() {
        return iterator;
    }

    public boolean hasNext() {
        return iterator.hasNext();
    }

    public Feature next() {
        return iterator.next();
    }

    public void remove() {
        iterator.remove();
    }
}


