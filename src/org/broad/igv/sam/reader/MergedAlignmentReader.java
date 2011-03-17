package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.io.IOException;
import java.util.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentQueryReader;

/**
 * Performs a logical merge of bam files.
 *
 * This implementation loads and closes iterators on each sam file sequentially to reduce the number of bam files
 * open at once.
 *
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
