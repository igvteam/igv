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

package org.broad.igv.sam.reader;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.sam.PicardAlignment;

/**
 *
 */
class SAMQueryIterator implements CloseableIterator<PicardAlignment> {

    String chr;
    int start;
    int end;
    boolean contained;
    /**
     * SAMRecord is 1-based, PicardAlignment is 0-based
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

    public PicardAlignment next() {
        SAMRecord ret = currentRecord;
        if (wrappedIterator.hasNext()) {
            currentRecord = wrappedIterator.next();
        } else {
            currentRecord = null;
        }
        return new PicardAlignment(ret);

    }

    public void remove() {
        //throw new UnsupportedOperationException("Not supported yet.");
    }
}
