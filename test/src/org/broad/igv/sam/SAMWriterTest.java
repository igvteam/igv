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

package org.broad.igv.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.MergedAlignmentReaderTest;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/05/04
 */
public class SAMWriterTest extends AbstractHeadlessTest {

    private String[] getTestPaths() {
        File indir = new File(TestUtils.DATA_DIR + "sam/");
        File[] infiles = indir.listFiles(new FileFilter() {
            @Override
            public boolean accept(File pathname) {
                String name = pathname.getAbsolutePath();
                return name.endsWith(".bam") || name.endsWith(".sam");
            }
        });

        List<String> inpaths = new ArrayList<String>(infiles.length);
        for (File f : infiles) {
            inpaths.add(f.getAbsolutePath());
        }

        return inpaths.toArray(new String[0]);
    }

    @Test
    public void testWriteRecordsFile() throws Exception {
        for (String inpath : getTestPaths()) {
            testWriteRecords(inpath, false);
        }
    }

    @Test
    public void testWriteRecordsStream() throws Exception {
        for (String inpath : getTestPaths()) {
            testWriteRecords(inpath, true);
        }
    }

    /**
     * Check that the alignments read from {@code outFile} are the same as those in {@code origAlignments}
     *
     * @param origAlignments
     * @param outFile
     * @param origPath       Used for error message only. Can be null
     * @throws java.io.IOException
     */
    public void checkRecordsMatch(List<PicardAlignment> origAlignments, File outFile, String origPath) throws IOException {
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
     * <p/>
     * If outStream is true, we use the stream writing methods of SAMWriter,
     * otherwise, file writing methods
     *
     * @throws Exception
     */
    public void testWriteRecords(String inpath, boolean outStream) throws Exception {
        String[] outpaths = new String[]{TestUtils.DATA_DIR + "out/tmp_sam.sam", TestUtils.DATA_DIR + "out/tmp_bam.bam"};
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
        private List<PicardAlignment> alignments;

        public SamHeaderIterator(String inpath) throws IOException {

            SAMReader reader = new SAMReader(inpath, false);
            this.header = reader.getFileHeader();
            Iterator<PicardAlignment> iter = reader.iterator();

            alignments = new ArrayList<PicardAlignment>();
            while (iter.hasNext()) {
                alignments.add(iter.next());
            }
        }

    }

    @Test
    public void testCopyBAMFile_01() throws Exception{
        String sequence = "chr1";
        int end = 300000000;
        int start = end / 5 - 1;
        String inpath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        ResourceLocator inlocator = new ResourceLocator(inpath);
        tstCopyBAMFile(inlocator, sequence, start, end);
    }

    @Test
    public void testCopyMergedBAM_01() throws Exception{
        String sequence = "chr1";
        int start = 151667156;
        int end = start + 10000;

        File listFile = new File(TestUtils.LARGE_DATA_DIR, "2largebams.bam.list");
        MergedAlignmentReaderTest.generateRepLargebamsList(listFile);
        ResourceLocator inlocator = new ResourceLocator(listFile.getAbsolutePath());
        tstCopyBAMFile(inlocator, sequence, start, end);
    }

    public void tstCopyBAMFile(ResourceLocator inlocator, String sequence, int start, int end) throws IOException{

        boolean createIndex = true;


        String outPath = TestUtils.TMP_OUTPUT_DIR + "tmpbam.bam";
        File outFile = new File(outPath);
        File indexFile = new File(outPath.replace(".bam", ".bai"));
        outFile.delete();
        indexFile.delete();
        outFile.deleteOnExit();
        indexFile.deleteOnExit();

        int writtenCount = SAMWriter.writeAlignmentFilePicard(inlocator, outPath, sequence, start, end);

        assertEquals("Index file existence unexpected: " + indexFile.getAbsolutePath(), createIndex, indexFile.exists());

        SAMFileReader writtenReader = new SAMFileReader(new File(outPath));
        writtenReader.setValidationStringency(ValidationStringency.SILENT);
        SAMRecordIterator iter = null;
        if(createIndex){
            iter = writtenReader.queryOverlapping(sequence, start + 1, end);
        }else{
            iter = writtenReader.iterator();
        }

        int readCount = 0;
        while (iter.hasNext()) {
            readCount++;
            iter.next();
        }

        System.out.println(readCount + " alignments read");
        assertTrue("No alignments read", readCount > 0);
        assertEquals("Read a different number of alignments than written", writtenCount, readCount);
    }
}
