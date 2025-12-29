package org.broad.igv.sam.reader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.logging.*;
import org.broad.igv.sam.EmptyAlignmentIterator;
import org.broad.igv.sam.SAMAlignment;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Sep 22, 2009
 * Time: 2:21:04 PM
 */
public class BAMReader implements AlignmentReader<SAMAlignment> {

    static Logger log = LogManager.getLogger(BAMReader.class);

    private SAMFileHeader header;
    private List<String> sequenceNames;
    private boolean indexed;
    private Map<String, Long> sequenceDictionary;
    private SamReaderPool readerPool;
    private int errorCount = 0;
    private String filename;

    public BAMReader(ResourceLocator locator, boolean requireIndex) throws IOException {

        indexed = requireIndex || locator.isHtsget();
        readerPool = new SamReaderPool(locator, requireIndex);
        filename = locator.getFileName();
        SamReader reader = readerPool.getReader();
        header = reader.getFileHeader();
        readerPool.freeReader(reader);
    }

    private SamReader getSamReader() throws IOException {
        return readerPool.getReader();
    }

    public void close() throws IOException {
        readerPool.close();
    }

    public synchronized SAMFileHeader getFileHeader() {
        if (header == null) {
            try {
                header = getSamReader().getFileHeader();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return header;
    }

    public boolean hasIndex() {
        return indexed;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getFileHeader());
    }

    public List<String> getSequenceNames() {
        if (sequenceNames == null) {
            loadSequenceDictionary();
        }
        return sequenceNames;
    }

    @Override
    public Map<String, Long> getSequenceDictionary() {
        if (sequenceDictionary == null) {
            loadSequenceDictionary();
        }
        return sequenceDictionary;
    }

    private void loadSequenceDictionary() {

        SAMFileHeader header = getFileHeader();
        sequenceNames = new ArrayList();
        sequenceDictionary = new HashMap<>();

        List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
        if (records.size() > 0) {
            for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                String chr = rec.getSequenceName();
                Long size = (long) rec.getSequenceLength();

                sequenceNames.add(chr);
                sequenceDictionary.put(chr, size);
            }
        }
    }


    public CloseableIterator<SAMAlignment> iterator() throws IOException {
        return new WrappedIterator(readerPool.getReaderIterator().iterator());
    }

    public CloseableIterator<SAMAlignment> query(String sequence, int start, int end, boolean contained) {
        if (sequenceDictionary != null && !sequenceDictionary.containsKey(sequence)) {
            return EMPTY_ITERATOR;
        } else {
            try {
                SamReader samReader = getSamReader();
                return new PicardIterator(samReader, sequence, start + 1, end, contained);
            } catch (Exception e) {
                if (errorCount == 0) {
                    MessageUtils.showMessage("Error querying alignments for: " + filename + "<br/>Error message: " + e.getMessage());
                }
                errorCount++;
                log.error("Error querying for sequence: " + sequence, e);
                return new EmptyAlignmentIterator();
            }
        }
    }

    static CloseableIterator<SAMAlignment> EMPTY_ITERATOR = new CloseableIterator<SAMAlignment>() {
        @Override
        public void close() {

        }

        @Override
        public boolean hasNext() {
            return false;
        }

        @Override
        public SAMAlignment next() {
            return null;
        }
    };


    class PicardIterator implements CloseableIterator<SAMAlignment> {

        private SamReader reader;    // Readers are 1-time use (for this iterator only)
        CloseableIterator<SAMRecord> iterator;

        public PicardIterator(SamReader samReader, String sequence, int start, int end, boolean contained) {
            this.reader = samReader;
            this.iterator = samReader.query(sequence, start, end, contained);
        }

        public void close() {
            iterator.close();
            readerPool.freeReader(this.reader);
        }

        public boolean hasNext() {
            return iterator.hasNext();
        }

        public SAMAlignment next() {
            return new SAMAlignment(iterator.next());
        }

        public void remove() {
            iterator.remove();
        }

    }
}
