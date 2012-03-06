package org.broad.igv.sam.reader;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SamAlignment;

/**
 *
 */
class SAMQueryIterator implements CloseableIterator<Alignment> {

    String chr;
    int start;
    int end;
    boolean contained;
    /**
     * SAMRecord is 1-based, SamAlignment is 0-based
     */
    SAMRecord currentRecord;
    CloseableIterator<SAMRecord> wrappedIterator;

    public SAMQueryIterator(CloseableIterator<SAMRecord> wrappedIterator) {
        this.chr = null;
        this.wrappedIterator = wrappedIterator;
        currentRecord = wrappedIterator.next();
    }

    /**
     *
     * @param sequence
     * @param start 0-based, inclusive-start
     * @param end 0-based, exclusive-end
     * @param contained
     * @param wrappedIterator
     */
    public SAMQueryIterator(String sequence, int start, int end, boolean contained,
                            CloseableIterator<SAMRecord> wrappedIterator) {
        this.chr = sequence;
        this.start = start;
        this.end = end;
        this.contained = contained;
        this.wrappedIterator = wrappedIterator;
        advanceToFirstRecord();
    }

    private void advanceToFirstRecord() {
        while (wrappedIterator.hasNext()) {
            currentRecord = wrappedIterator.next();
            if (!currentRecord.getReferenceName().equals(chr)) {
                break;
            //currentRecord is 1-based, end-inclusive.
            //start/end are 0-based, end-exclusive
            } else if ((contained && currentRecord.getAlignmentStart()-1 >= start) ||
                    (!contained && currentRecord.getAlignmentEnd()-1 >= start)) {
                break;
            }
        }
    }

    public void close() {
        wrappedIterator.close();
    }

    public boolean hasNext() {
        if (chr == null && currentRecord != null) {
            return true;
        }
        if (currentRecord == null || (chr != null && !chr.equals(currentRecord.getReferenceName()))) {
            return false;
        } else {
            return contained ? currentRecord.getAlignmentEnd() <= end
                    : currentRecord.getAlignmentStart() <= end;
        }
    }

    public SamAlignment next() {
        SAMRecord ret = currentRecord;
        if (wrappedIterator.hasNext()) {
            currentRecord = wrappedIterator.next();
        } else {
            currentRecord = null;
        }
        return new SamAlignment(ret);

    }

    public void remove() {
        //throw new UnsupportedOperationException("Not supported yet.");
    }
}
