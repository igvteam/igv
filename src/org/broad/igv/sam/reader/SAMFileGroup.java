/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

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

