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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * @author jrobinso
 */
public class SAMFileGroup implements Iterable<SAMRecord> {

    // Map of chromosome -> file
    File[] samFiles;

    public SAMFileGroup(File[] samFiles) {
        this.samFiles = samFiles;
    }

    public void close() throws IOException {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public SAMFileHeader getHeader() {
        if (samFiles.length == 0) {
            return null;
        }
        SAMFileReader reader = new SAMFileReader(samFiles[0]);
        SAMFileHeader header = reader.getFileHeader();
        reader.close();
        ;
        return header;
    }

    public Iterator<SAMRecord> iterator() {
        return new SamGroupIterator();
    }

    public class SamGroupIterator implements CloseableIterator<SAMRecord> {

        int idx;
        SAMFileReader currentReader;
        CloseableIterator<SAMRecord> currentIter;
        SAMRecord next;

        SamGroupIterator() {
            idx = 0;
            advance();
        }

        public void close() {
            throw new UnsupportedOperationException("Not supported yet.");
        }

        public boolean hasNext() {
            return next != null;
        }

        public SAMRecord next() {
            SAMRecord ret = next;
            advance();
            return ret;
        }

        private void advance() {
            if (currentIter == null || !currentIter.hasNext()) {
                if (currentIter != null) {
                    currentIter.close();
                    currentIter = null;
                }
                if (currentReader != null) {
                    currentReader.close();
                }
                if (idx < samFiles.length) {
                    currentReader = new SAMFileReader(samFiles[idx]);
                    currentIter = currentReader.iterator();
                    idx++;
                }
            }
            next = currentIter == null ? null : currentIter.next();
        }

        public void remove() {
            throw new UnsupportedOperationException("Not supported yet.");
        }
    }
}

