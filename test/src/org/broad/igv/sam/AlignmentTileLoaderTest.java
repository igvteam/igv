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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class AlignmentTileLoaderTest extends AbstractHeadlessTest {



    /**
     * Test that sampling keeps pairs together.
     *
     * @throws Exception
     */
    @Test
    public void testKeepPairsDownsample_02() throws Exception {
        String path = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";

        String sequence = "chr1";
        int start = 151766945;
        int end = 151791197;
        int maxDepth = 5;

        AlignmentTileLoader.AlignmentTile tile= tstKeepPairsDownsample(path, sequence, start, end, maxDepth);
        assertTrue(tile.getDownsampledIntervals().size() > 0);
    }

    /**
     * Test that sampling keeps pairs together.
     *
     * @throws Exception
     */
    @Test
    public void testNoDownsample() throws Exception {
        String path = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";

        String sequence = "chr1";
        int start = 151766945;
        int end = 151791197;
        int maxDepth = 500;

        AlignmentTileLoader.AlignmentTile tile= tstKeepPairsDownsample(path, sequence, start, end, maxDepth);
        assertEquals(0, tile.getDownsampledIntervals().size());


    }

    private AlignmentTileLoader.AlignmentTile tstKeepPairsDownsample(String path, String sequence, int start, int end, int maxDepth) throws Exception{


        String oldMaxVis = PreferenceManager.getInstance().get(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "" + (end - start));

        int actMaxDepth = 100;
        if(maxDepth > 0){
            actMaxDepth = maxDepth;
        }

        try {
            ResourceLocator loc = new ResourceLocator(path);
            AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
            AlignmentTileLoader loader = new AlignmentTileLoader(reader);

            AlignmentDataManager.DownsampleOptions downsampleOptions = new AlignmentDataManager.DownsampleOptions(true, 50, actMaxDepth);

            AlignmentTileLoader.AlignmentTile tile = loader.loadTile(sequence, start, end, null, downsampleOptions, null, null, null);
            List<Alignment> alignments = tile.getAlignments();
            int count = 0;
            Map<String, Integer> pairedReads = new HashMap<String, Integer>();
            for(Alignment al: alignments) {
                assertNotNull(al);
                count++;

                //Only look at proper pairs, which are a subset.
                //Our system should keep all things with the same read name, but we don't know
                //how many chimeric/secondary alignments there might be
                if (al.isProperPair()) {
                    //Mate may not be part of the query.
                    //Make sure it's within bounds
                    int mateStart = al.getMate().getStart();
                    //All we require is some overlap
                    boolean overlap = (mateStart + al.getReadSequence().length()) >= start && mateStart < end;
                    overlap &= al.getMate().getChr().equals(al.getChr());
                    if (overlap) {
                        Integer rdCnt = pairedReads.get(al.getReadName());
                        rdCnt = rdCnt != null ? rdCnt + 1 : 1;
                        pairedReads.put(al.getReadName(), rdCnt);
                    }
                }
            }

            assertTrue("No alignments loaded", count > 0);

            int countmissing = 0;
            for (String readName : pairedReads.keySet()) {
                int val = pairedReads.get(readName);
                countmissing += 2 == val ? 0 : 1;
                if (val != 2) {
                    System.out.println("Read " + readName + " has val " + val);
                }
            }

            System.out.println("Number of paired reads: " + pairedReads.size());
            assertTrue("No pairs in test data set", pairedReads.size() > 0);
            assertEquals("Missing " + countmissing + " out of " + pairedReads.size() + " pairs", 0, countmissing);
            return tile;
        } catch (Exception e) {
            throw e;
        }finally{
            PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, oldMaxVis);
        }

    }

}