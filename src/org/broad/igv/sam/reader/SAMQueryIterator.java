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
    SAMRecord currentRecord;
    CloseableIterator<SAMRecord> wrappedIterator;

    public SAMQueryIterator(CloseableIterator<SAMRecord> wrappedIterator) {
        this.chr = null;
        this.wrappedIterator = wrappedIterator;
        currentRecord = wrappedIterator.next();
    }

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
            } else if ((contained && currentRecord.getAlignmentStart() >= start) ||
                    (!contained && currentRecord.getAlignmentEnd() >= start)) {
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
