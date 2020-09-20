package org.broad.igv.sam.reader;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.zip.InflaterFactory;
import org.apache.log4j.Logger;
import org.broad.igv.sam.EmptyAlignmentIterator;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.*;

import static org.broad.igv.sam.reader.BAMReader.EMPTY_ITERATOR;

public class HtsgetBAMReader implements AlignmentReader<PicardAlignment> {

    static Logger log = Logger.getLogger(HtsgetBAMReader.class);

    private final ResourceLocator locator;

    HtsgetBAMFileReader reader;
    SAMFileHeader header;
    List<String> sequenceNames;
    private boolean indexed = false; // False until proven otherwise
    private Map<String, Long> sequenceDictionary;

    public HtsgetBAMReader(ResourceLocator locator, boolean requireIndex) throws IOException, URISyntaxException {
        this.locator = locator;
        reader = getHtsgetReader(locator, requireIndex);
        header = reader.getFileHeader();
    }

    private HtsgetBAMFileReader getHtsgetReader(ResourceLocator locator, boolean requireIndex) throws URISyntaxException, IOException {
        return new HtsgetBAMFileReader(
                new URI(this.locator.getURLPath()),
                false,
                ValidationStringency.SILENT,
                DefaultSAMRecordFactory.getInstance(),
                false
        );
    }

    @Override
    public void close() throws IOException {
        if (reader != null) {
            reader.close();
        }
    }

    @Override
    public List<String> getSequenceNames() throws IOException {
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
                Long size = new Long(rec.getSequenceLength());

                sequenceNames.add(chr);
                sequenceDictionary.put(chr, size);
            }
        }

    }

    @Override
    public SAMFileHeader getFileHeader() {
        if (header == null) {
            header = reader.getFileHeader();
        }
        return header;
    }

    @Override
    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getFileHeader());
    }

    @Override
    public CloseableIterator<PicardAlignment> iterator() {
        return new WrappedIterator(reader.getIterator());
    }

    @Override
    public CloseableIterator<PicardAlignment> query(String sequence, int start, int end, boolean contained) throws IOException {
        if (sequenceDictionary != null && !sequenceDictionary.containsKey(sequence)) {
            return EMPTY_ITERATOR;
        } else {
            CloseableIterator<SAMRecord> iter = null;
            try {
                iter = reader.query(sequence, start + 1, end, contained);
            } catch (IllegalArgumentException e) {
                log.error("Error querying for sequence: " + sequence, e);
                return new EmptyAlignmentIterator();
            }
            return new BAMReader.ListIterator(iter);
        }
    }

    @Override
    public boolean hasIndex() {
        return indexed;
    }
}
