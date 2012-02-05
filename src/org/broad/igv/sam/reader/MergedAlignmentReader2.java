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
import net.sf.samtools.util.CloseableIterator;

import java.io.IOException;
import java.util.*;

import org.broad.igv.sam.Alignment;

/**
 * Performs a logical merge of bam files.
 *
 * This implementation loads and closes iterators on each sam file sequentially to reduce the number of bam files
 * open at once.
 *
 * User: jrobinso
 * Date: Apr 25, 2010
 */
public class MergedAlignmentReader2 implements AlignmentReader {

    Collection<AlignmentReader> readers;

    public MergedAlignmentReader2(Collection<AlignmentReader> readers) {
        this.readers = readers;
    }

    public CloseableIterator<Alignment> iterator() {
        return new MergedFileIterator();
    }

    public CloseableIterator<Alignment> query(String chr, int start, int end, boolean contained) throws IOException {
        return new MergedFileIterator(chr, start, end, contained);
    }

    public void close() throws IOException {
        for (AlignmentReader reader : readers) {
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
            for (AlignmentReader reader : readers) {
                CloseableIterator<Alignment> iter = reader.iterator();
                if (iter.hasNext()) {
                    iterators.add(new RecordIterWrapper(iter));
                }
                iter.close();
            }
        }

        public MergedFileIterator(String chr, int start, int end, boolean contained) throws IOException {
            iterators = new PriorityQueue(readers.size(), new FooComparator());
            for (AlignmentReader reader : readers) {
                CloseableIterator<Alignment> iter = reader.query(chr, start, end, contained);
                if (iter.hasNext()) {
                    ArrayList<Alignment> records = new ArrayList();
                    while (iter.hasNext()) {
                        records.add(iter.next());
                    }
                    iter.close();
                    iterators.add(new RecordIterWrapper(records.iterator()));
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
            // no-op
        }

        class RecordIterWrapper {

            Alignment nextRecord;
            Iterator<Alignment> iterator;

            RecordIterWrapper(Iterator<Alignment> iter) {
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
        }

        class FooComparator implements Comparator<RecordIterWrapper> {

            public int compare(RecordIterWrapper wrapper, RecordIterWrapper foo1) {
                return wrapper.nextRecord.getAlignmentStart() - foo1.nextRecord.getAlignmentStart();
            }
        }
    }
}
