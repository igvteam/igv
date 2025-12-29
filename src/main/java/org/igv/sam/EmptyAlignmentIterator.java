/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.sam;

import htsjdk.samtools.util.CloseableIterator;

public class EmptyAlignmentIterator implements CloseableIterator<SAMAlignment> {

    static EmptyAlignmentIterator instance = new EmptyAlignmentIterator();

    public static EmptyAlignmentIterator getInstance() {
        return instance;
    }

    public EmptyAlignmentIterator() {
    }

    public void close() {
        // ignore
    }

    public boolean hasNext() {
        return false;
    }

    public SAMAlignment next() {
        return null;
    }

    public void remove() {
        // ignore
    }
}
