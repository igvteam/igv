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
import org.broad.igv.Globals;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.BAMRemoteQueryReader;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class BAMCachingQueryReaderTest {

    //String testFile = "/xchip/igv/testfiles/sam/NA12878.chrom1.SRP000032.2009_02.bam";
    String testFile = "/xchip/igv/data/public/BodyMap/hg18/50bp/FCA/s_1_1_sequence.bam";
    String server = "http://www.broadinstitute.org/webservices/igv";
    String sequence = "chr1";
    int start = 44780145 - 100000;
    int end = 44789983;
    //String testFile = "/Volumes/igv/testfiles/Sam/303KY.8.paired.bam";
    //String server = null;
    //String sequence = "chr1";    //int start = 58892867;
    //int end = 58960732;

    public BAMCachingQueryReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Globals.setHeadless(true);
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
     * Test of getHeader method, of class BAMCachingQueryReader. The test compares
     * the results of BAMCachingQueryReader with SamRemoteQueryReader, which
     * is assumed to be correct.
     */
    @Test
    public void testGetHeader() throws IOException {
        ResourceLocator loc = new ResourceLocator(server, testFile);
        AlignmentQueryReader reader = new BAMRemoteQueryReader(loc);
        SAMFileHeader expectedHeader = reader.getHeader();
        reader.close();


        CachingQueryReader cachingReader = new CachingQueryReader(new BAMRemoteQueryReader(loc));
        SAMFileHeader header = cachingReader.getHeader();
        cachingReader.close();

        assertTrue(header.equals(expectedHeader));


    }

    /**
     * Test of query method, of class BAMCachingQueryReader.  The test compares
     * the results of CachingQueryReader with BAMRemoteQueryReader, which
     * is assumed to be correct.
     */
    @Test
    public void testQuery() throws IOException {
        boolean contained = false;

        ResourceLocator loc = new ResourceLocator(server, testFile);
        AlignmentQueryReader reader = new BAMRemoteQueryReader(loc);
        CloseableIterator<Alignment> iter = reader.query(sequence, start, end, contained);
        //TODO the results may be returned in a different order. Not sure if that's a bug or not
        Map<String, Alignment> expectedResult = new HashMap();
        while (iter.hasNext()) {
            Alignment rec = iter.next();
            expectedResult.put(rec.getReadName(), rec);
        }
        reader.close();

        CachingQueryReader cachingReader = new CachingQueryReader(new BAMRemoteQueryReader(loc));
        CloseableIterator<Alignment> cachingIter = cachingReader.query(sequence, start, end, new ArrayList(), 100, null);
        List<Alignment> result = new ArrayList();

        while (cachingIter.hasNext()) {
            result.add(cachingIter.next());
        }
        cachingReader.close();

        assertTrue(expectedResult.size() > 0);
        assertEquals(expectedResult.size(), result.size());
        for (int i = 0; i < result.size(); i++) {
            Alignment res = result.get(i);
            assertTrue(expectedResult.containsKey(res.getReadName()));
            Alignment exp = expectedResult.get(res.getReadName());
            assertEquals(exp.getAlignmentStart(), res.getAlignmentStart());
            assertEquals(exp.getAlignmentEnd(), res.getAlignmentEnd());
        }
    }
}