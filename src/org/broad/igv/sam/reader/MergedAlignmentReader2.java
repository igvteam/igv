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
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.util.CloseableIterator;

import java.io.IOException;
import java.util.*;

import org.broad.igv.sam.Alignment;

/**
 * Performs a logical merge of bam files.
 * <p/>
 * This implementation loads and closes iterators on each sam file sequentially to reduce the number of bam files
 * open at once.
 * <p/>
 * User: jrobinso
 * Date: Apr 25, 2010
 */
public class MergedAlignmentReader2 implements AlignmentReader {

    List<AlignmentReader> readers;
    List<String> sequenceNames;
    Map<String, Integer> sequenceNameIndex;

    public MergedAlignmentReader2(List<AlignmentReader> readers) {
        this.readers = readers;
        loadSequenceNames();
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


    public List<String> getSequenceNames() {
        return sequenceNames;
    }

    public Set<String> getPlatforms() {
        Set<String> platforms = new HashSet<String>();
        for (AlignmentReader reader : readers) {
            Set<String> plf = reader.getPlatforms();
            if(plf != null){
                platforms.addAll(plf);
            }
        }
        return platforms;
    }


    /**
     * Return the merged list of all sequence names, maintaining order.
     *
     * @return
     */
    public void loadSequenceNames() {
        // Use a set for quick comparison
        LinkedHashSet<String> names = new LinkedHashSet<String>(50);
        for (AlignmentReader reader : readers) {
            names.addAll(reader.getSequenceNames());
        }
        sequenceNames = new ArrayList<String>(names);
        sequenceNameIndex = new HashMap<String, Integer>(sequenceNames.size());
        for (int i = 0; i < sequenceNames.size(); i++) {
            sequenceNameIndex.put(sequenceNames.get(i), i);
        }
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
            iterators = new PriorityQueue(readers.size(), new AlignmentStartComparator());
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

        class AlignmentStartComparator implements Comparator<RecordIterWrapper> {

            public int compare(RecordIterWrapper wrapper1, RecordIterWrapper wrapper2) {
                Alignment a1 = wrapper1.nextRecord;
                Alignment a2 = wrapper2.nextRecord;

                int idx1 = sequenceNameIndex.get(a1.getChr());
                int idx2 = sequenceNameIndex.get(a2.getChr());
                if (idx1 > idx2) {
                    return 1;
                } else if (idx1 < idx2) {
                    return -1;
                } else {
                    return a1.getAlignmentStart() - a2.getAlignmentStart();
                }
            }
        }
    }
}