/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
import org.broad.igv.sam.Alignment;

import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Set;


/**
 * Performs a logical merge of bam files
 * User: jrobinso
 * Date: Apr 25, 2010
 */
public class MergedAlignmentReader implements AlignmentQueryReader {

    Collection<AlignmentQueryReader> readers;

    public MergedAlignmentReader(Collection<AlignmentQueryReader> readers) {
        this.readers = readers;
    }

    public CloseableIterator<Alignment> iterator() {
        return new MergedFileIterator();
    }

    public CloseableIterator<Alignment> query(String chr, int start, int end, boolean contained) throws IOException {
        return new MergedFileIterator(chr, start, end, contained);
    }

    public void close() throws IOException {
        for (AlignmentQueryReader reader : readers) {
            reader.close();
        }
    }


    /**
     * Return the header for this file group.
     * // TODO -- we are assuming that the headers are consistent for all files in this group, we should really check
     *
     * @return
     */
    public SAMFileHeader getHeader() throws IOException {
        return readers.iterator().next().getHeader();
    }

    public Set<String> getSequenceNames() {
        return readers.iterator().next().getSequenceNames();
    }

    public boolean hasIndex() {
        return readers.iterator().next().hasIndex();
    }


    public class MergedFileIterator implements CloseableIterator<Alignment> {

        PriorityQueue<RecordIterWrapper> iterators;

        public MergedFileIterator() {
            iterators = new PriorityQueue(readers.size());
            for (AlignmentQueryReader reader : readers) {
                CloseableIterator<Alignment> iter = reader.iterator();
                if (iter.hasNext()) {
                    iterators.add(new RecordIterWrapper(iter));
                }
                iter.close();
            }
        }

        public MergedFileIterator(String chr, int start, int end, boolean contained) throws IOException {
            iterators = new PriorityQueue(readers.size(), new FooComparator());
            for (AlignmentQueryReader reader : readers) {
                CloseableIterator<Alignment> iter = reader.query(chr, start, end, contained);
                if (iter.hasNext()) {
                    iterators.add(new RecordIterWrapper(iter));
                }
            }
        }

        public boolean hasNext() {
            return iterators.size() > 0;
        }

        public Alignment next() {
            RecordIterWrapper wrapper = iterators.poll();
            Alignment next = wrapper.advance();
            if (wrapper.hasNext()) {
                iterators.add(wrapper);
            } else {
                //   wrapper.close();
            }
            return next;
        }

        public void remove() {
            throw new UnsupportedOperationException("Remove not implemented");
        }

        public void close() {
            for (RecordIterWrapper wrapper : iterators) {
                wrapper.close();
            }
        }

        class RecordIterWrapper {

            Alignment nextRecord;
            CloseableIterator<Alignment> iterator;

            RecordIterWrapper(CloseableIterator<Alignment> iter) {
                this.iterator = iter;
                nextRecord = (iterator.hasNext() ? iterator.next() : null);
            }

            Alignment advance() {
                Alignment tmp = nextRecord;
                nextRecord = (iterator.hasNext() ? iterator.next() : null);
                return tmp;

            }

            boolean hasNext() {
                return nextRecord != null;
            }

            void close() {
                if (iterator != null) {
                    iterator.close();
                    iterator = null;
                }
            }
        }

        class FooComparator implements Comparator<RecordIterWrapper> {

            public int compare(RecordIterWrapper wrapper, RecordIterWrapper foo1) {
                return wrapper.nextRecord.getAlignmentStart() - foo1.nextRecord.getAlignmentStart();
            }
        }
    }
}
