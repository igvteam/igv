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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.ReadGroupFilter;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class CachingQueryReaderTest {

    String testFile = "http://www.broadinstitute.org/igvdata/BodyMap/hg18/50bp/FCA/s_1_1_sequence.bam";
    String sequence = "chr1";
    int start = 44680145;
    int end = 44789983;
    private boolean contained = false;
    ;

    public CachingQueryReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        TestUtils.setUpHeadless();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of getHeader method, of class CachingQueryReader. The test compares
     * the results of CachingQueryReader with a non-caching reader which
     * is assumed to be correct.
     */
    @Test
    public void testGetHeader() throws IOException {

        ResourceLocator loc = new ResourceLocator(testFile);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        SAMFileHeader expectedHeader = reader.getHeader();
        reader.close();

        reader = AlignmentReaderFactory.getReader(loc);
        CachingQueryReader cachingReader = new CachingQueryReader(reader);
        SAMFileHeader header = cachingReader.getHeader();
        cachingReader.close();

        assertTrue(header.equals(expectedHeader));
    }

    @Test
    public void testQuery() throws IOException {
        tstQuery(testFile, sequence, start, end, contained, Integer.MAX_VALUE / 1000);
    }

    /**
     * Test of query method, of class CachingQueryReader.  The test compares
     * the results of CachingQueryReader non-caching reader which
     * is assumed to be correct.
     *
     * Note that SAMFileReader (which is the non-caching reader) is 1-based
     * and inclusive-end. CachingQueryReader is 0-based and exclusive end.
     */
    public void tstQuery(String testFile, String sequence, int start, int end, boolean contained, int maxDepth) throws IOException {

        ResourceLocator loc = new ResourceLocator(testFile);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        CloseableIterator<Alignment> iter = reader.query(sequence, start+1, end, contained);

        Map<String, Alignment> expectedResult = new HashMap();
        while (iter.hasNext()) {
            Alignment rec = iter.next();

            // the following filters are applied in the Caching reader, so we need to apply them here.
            boolean filterFailedReads = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_FILTER_FAILED_READS);
            ReadGroupFilter filter = ReadGroupFilter.getFilter();
            boolean showDuplicates = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_DUPLICATES);
            int qualityThreshold = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_QUALITY_THRESHOLD);
            if (!rec.isMapped() || (!showDuplicates && rec.isDuplicate()) ||
                    (filterFailedReads && rec.isVendorFailedRead()) ||
                    rec.getMappingQuality() < qualityThreshold ||
                    (filter != null && filter.filterAlignment(rec))) {
                continue;
            }

            expectedResult.put(rec.getReadName(), rec);
            if(contained){
                assertTrue(rec.getStart() >= start);
            }else{
                //All we require is some overlap
                boolean overlap = rec.getStart() >= start && rec.getStart() < end;
                overlap |= (rec.getEnd() >= start) && (rec.getStart() < start);
                assertTrue(overlap);
            }
            assertEquals(sequence, rec.getChr());
        }
        reader.close();

        reader = AlignmentReaderFactory.getReader(loc);
        CachingQueryReader cachingReader = new CachingQueryReader(reader);
        CloseableIterator<Alignment> cachingIter = cachingReader.query(sequence, start, end, new ArrayList(),
                new ArrayList(), maxDepth, null, null);
        List<Alignment> result = new ArrayList();

        while (cachingIter.hasNext()) {
            result.add(cachingIter.next());
        }
        cachingReader.close();


        assertTrue(expectedResult.size() > 0);
        assertEquals(expectedResult.size(), result.size());
        for (int i = 0; i < result.size(); i++) {
            Alignment rec = result.get(i);

            if(contained){
                assertTrue(rec.getStart() >= start);
            }else{
                //All we require is some overlap
                boolean overlap = rec.getStart() >= start && rec.getStart() <= end;
                overlap |= start >= rec.getStart() && start <= rec.getEnd();
                assertTrue(overlap);
            }
            assertEquals(sequence, rec.getChr());

            assertTrue(expectedResult.containsKey(rec.getReadName()));
            Alignment exp = expectedResult.get(rec.getReadName());
            assertEquals(exp.getAlignmentStart(), rec.getAlignmentStart());
            assertEquals(exp.getAlignmentEnd(), rec.getAlignmentEnd());
        }
    }

    //@Ignore
    @Test
    public void testQueryLargeFile() throws Exception{
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "5");
        String path = TestUtils.LARGE_DATA_DIR + "/ABCD_igvSample.bam";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        CachingQueryReader cachingReader = new CachingQueryReader(reader);

        //Edge location
        String sequence = "chr12";
        int start = 56815621-1;
        int end = start+1;
        int expSize = 1066;

        tstSize(cachingReader, sequence,  start, end, expSize * 5, expSize);
        tstQuery(path, sequence,  start, end, false, 10000);

        //Edge location
        sequence = "chr12";
        start = 56815644-1;
        end = start;
        expSize = 271;

        //tstSize(cachingReader, sequence,  start, end, expSize * 5, expSize);
        tstQuery(path, sequence,  start, end, false, 10000);

        //Center location
        sequence = "chr12";
        start = 56815674;
        end = start+1;

        expSize = 3288;

        //tstSize(cachingReader, sequence,  start, end, expSize * 5, expSize);
        tstQuery(path, sequence,  start, end, false, 10000);


    }

    @Test
    public void testQueryPiledUp() throws Exception{
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "5");
        String path = TestUtils.DATA_DIR + "/aligned/pileup.sorted.aligned";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        CachingQueryReader cachingReader = new CachingQueryReader(reader);

        //Edge location
        String sequence = "chr1";
        int start = 141;
        int end = start+1;
        int expSize = 40;

        tstSize(cachingReader, sequence,  start, end, expSize * 5, expSize);

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        cachingReader = new CachingQueryReader(reader);

        tstSize(cachingReader, sequence,  start, end, expSize * 100, expSize);
        tstQuery(path, sequence,  start, end, false, 10000);

        //Center, deep coverage region
        sequence = "chr1";
        start = 429;
        end = start+1;
        int coverageLim = 1000;
        expSize = 1408;

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        cachingReader = new CachingQueryReader(reader);


        //tstSize(cachingReader, sequence,  start, end, coverageLim, expSize);

        coverageLim = 10000;
        expSize = 1408;

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        cachingReader = new CachingQueryReader(reader);


        //tstSize(cachingReader, sequence,  start, end, coverageLim, expSize);

        //This doesn't work on .aligned files, the query returns improper results
        tstQuery(path, sequence,  start, end, false, coverageLim);

    }
    
    public List<Alignment> tstSize(CachingQueryReader cachingReader, String sequence, int start, int end, int maxDepth, int expSize){
        CloseableIterator<Alignment> cachingIter = cachingReader.query(sequence, start, end, new ArrayList(),
                new ArrayList(), maxDepth, null, null);
        List<Alignment> result = new ArrayList();

        while (cachingIter.hasNext()) {
            result.add(cachingIter.next());
        }
        
        assertEquals(expSize, result.size());
        return result;
    }

    /**
     * The main purpose of this test is to see if we get a
     * heap space error.
     *
     * @throws Exception
     */
    @Test
    public void testQueryLargeFile2() throws Exception{
        String path = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.454.bam";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        CachingQueryReader cachingReader = new CachingQueryReader(reader);

        String sequence = "MT";
        int start = 1000;
        int end = 3000;
        int maxDepth = 1000;

        CloseableIterator<Alignment> iter = cachingReader.query(sequence,  start, end, new ArrayList(),
                new ArrayList(), maxDepth, null, null);
        int count=0;
        while(iter.hasNext()){
            assertNotNull(iter.next());
            count++;
        }
        
        assertTrue(count > 0);


    }

}