package org.igv.alignment;

import htsjdk.samtools.*;
import org.igv.AbstractHeadlessTest;
import org.igv.alignment.reader.AlignmentReader;
import org.igv.alignment.reader.AlignmentReaderFactory;
import org.igv.alignment.reader.SAMReader;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012/05/04
 */
public class SAMWriterTest extends AbstractHeadlessTest {

    final static String ifile = TestUtils.DATA_DIR + "/sam/NA12878.muc1.test.sam";

    @Test
    public void testWriteRecordsFile() throws Exception {
        testWriteRecords(ifile, false);

    }

    @Test
    public void testWriteRecordsStream() throws Exception {
        testWriteRecords(ifile, true);
    }

    /**
     * Check that the alignments read from {@code outFile} are the same as those in {@code origAlignments}
     *
     * @param origAlignments
     * @param outFile
     * @param origPath       Used for error message only. Can be null
     * @throws java.io.IOException
     */
    public void checkRecordsMatch(List<SAMAlignment> origAlignments, File outFile, String origPath) throws IOException {
        //Read back in, check equality
        AlignmentReader outputReader = AlignmentReaderFactory.getReader(outFile.getAbsolutePath(), false);
        Iterator<Alignment> outputIter = outputReader.iterator();
        int index = 0;
        while (outputIter.hasNext()) {
            Alignment outputAl = outputIter.next();
            Alignment origAl = origAlignments.get(index);

            String errmsg = "Cigar strings not equal at " + index + " in " + origAl.getReadName();
            if (origPath != null) {
                errmsg += " in file " + origPath;
            }
            assertEquals(errmsg, origAl.getCigarString(), outputAl.getCigarString());
            index++;
        }
        assertEquals("Incorrect number of alignments", origAlignments.size(), index);
    }

    /**
     * Test our ability to write SAM Records out.
     * We check both bam and sam output format
     * <p>
     * If outStream is true, we use the stream writing methods of SAMWriter,
     * otherwise, file writing methods
     *
     * @throws Exception
     */
    public void testWriteRecords(String inpath, boolean outStream) throws Exception {
        String[] outpaths = new String[]{TestUtils.TMP_OUTPUT_DIR + "tmp_sam.sam", TestUtils.TMP_OUTPUT_DIR + "tmp_bam.bam"};
        for (String outpath : outpaths) {
            SamHeaderIterator shi = new SamHeaderIterator(inpath);

            File outFile = new File(outpath);
            outFile.delete();
            outFile.deleteOnExit();


            SAMWriter writer = new SAMWriter(shi.header);

            boolean bam = outpath.endsWith(".bam");

            if (!outStream) {
                writer.writeToFile(outFile, shi.alignments.iterator(), false);
            } else {
                OutputStream outputStream = new BufferedOutputStream(new FileOutputStream(outFile));
                writer.writeToStream(outputStream, shi.alignments.iterator(), bam);
            }
            checkRecordsMatch(shi.alignments, outFile, inpath);

        }
    }

    private static class SamHeaderIterator {
        private SAMFileHeader header;
        private List<SAMAlignment> alignments;

        public SamHeaderIterator(String inpath) throws IOException {

            SAMReader reader = new SAMReader(inpath, false);
            this.header = reader.getFileHeader();
            Iterator<SAMAlignment> iter = reader.iterator();

            alignments = new ArrayList<SAMAlignment>();
            while (iter.hasNext()) {
                alignments.add(iter.next());
            }
        }

    }
}
