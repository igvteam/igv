package org.igv.feature;

import htsjdk.tribble.CloseableTribbleIterator;

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


