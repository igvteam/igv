/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;
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
    public void checkRecordsMatch(List<SamAlignment> origAlignments, File outFile, String origPath) throws IOException {
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
        private List<SamAlignment> alignments;

        public SamHeaderIterator(String inpath) throws IOException {

            SAMReader reader = new SAMReader(inpath, false);
            this.header = reader.getFileHeader();
            Iterator<SamAlignment> iter = reader.iterator();

            alignments = new ArrayList<SamAlignment>();
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
        writtenReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
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
