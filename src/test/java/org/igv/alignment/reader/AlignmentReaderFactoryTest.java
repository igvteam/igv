package org.igv.alignment.reader;

import htsjdk.samtools.BamIndexType;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.igv.AbstractHeadlessTest;
import org.igv.alignment.Alignment;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class AlignmentReaderFactoryTest extends AbstractHeadlessTest {

    @Test
    public void bgzippedSamWithCsiIndexIsQueriedThroughIt() throws Exception {
        File sam = writeIndexedSam("csi.sam.gz", BamIndexType.CSI);
        assertTrue(new File(sam.getPath() + ".csi").exists());

        AlignmentReader reader = AlignmentReaderFactory.getReader(sam.getPath(), true);

        assertTrue(reader.hasIndex());
        assertEquals(List.of("r2"), readNames(reader.query("chr1", 1500, 2500, false)));
    }

    @Test
    public void bgzippedSamWithBaiIndexIsQueriedThroughIt() throws Exception {
        File sam = writeIndexedSam("bai.sam.gz", BamIndexType.BAI);
        assertTrue(new File(sam.getPath() + ".bai").exists());

        AlignmentReader reader = AlignmentReaderFactory.getReader(sam.getPath(), true);

        assertTrue(reader.hasIndex());
        assertEquals(List.of("r2"), readNames(reader.query("chr1", 1500, 2500, false)));
    }

    @Test
    public void bgzippedSamIsReadThroughWithoutAnIndex() throws Exception {
        File sam = writeIndexedSam("unindexed.sam.gz", BamIndexType.CSI);
        new File(sam.getPath() + ".csi").delete();

        AlignmentReader reader = AlignmentReaderFactory.getReader(sam.getPath(), false);

        assertEquals(List.of("r1", "r2", "r3"), readNames(reader.iterator()));
    }

    /**
     * Writes a bgzipped SAM with 10-base reads r1, r2 and r3 on chr1 at 1000, 2000 and 3000, indexed while writing.
     */
    private static File writeIndexedSam(String name, BamIndexType indexType) {
        File dir = new File(TestUtils.TMP_OUTPUT_DIR, "alignment_reader_factory");
        dir.mkdirs();
        File sam = new File(dir, name);

        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 100000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true).setSamIndexType(indexType);
        try (SAMFileWriter writer = factory.makeWriter(header, true, sam.toPath(), null)) {
            for (int i = 1; i <= 3; i++) {
                SAMRecord record = new SAMRecord(header);
                record.setReadName("r" + i);
                record.setReferenceName("chr1");
                record.setAlignmentStart(i * 1000);
                record.setMappingQuality(60);
                record.setCigarString("10M");
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                writer.addAlignment(record);
            }
        }
        return sam;
    }

    private static List<String> readNames(CloseableIterator<? extends Alignment> alignments) {
        List<String> names = new ArrayList<>();
        try (alignments) {
            alignments.forEachRemaining(a -> names.add(a.getReadName()));
        }
        return names;
    }
}
