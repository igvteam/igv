package org.igv.sam.reader;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.igv.sam.SAMAlignment;

/**
 * @author jrobinso
 * @since Sep 22, 2009
 */
public class WrappedIterator implements CloseableIterator<SAMAlignment> {

    CloseableIterator<SAMRecord> iter;

    public WrappedIterator(CloseableIterator<SAMRecord> iter) {
        this.iter = iter;
    }

    public void close() {
        iter.close();
    }

    public boolean hasNext() {
        return iter.hasNext();
    }

    public SAMAlignment next() {
        return new SAMAlignment(iter.next());
    }

    public void remove() {
        iter.remove();
    }
}
