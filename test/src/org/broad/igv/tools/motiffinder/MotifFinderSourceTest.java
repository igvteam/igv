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

package org.broad.igv.tools.motiffinder;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.util.Iterator;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2013-Jan-22
 */
public class MotifFinderSourceTest extends AbstractHeadlessTest{


    String shortSeq = "GACTCTGACTGGACTCTGATCAG";

    @Test
    public void testBasicSearch() throws Exception{
        String motif = "TCTG";
        int posStart = 19389;
        Iterator<Feature> matchIter = MotifFinderSource.search(null, motif, posStart, shortSeq.getBytes());
        assertTrue(matchIter.hasNext());
        Feature feat = matchIter.next();

        assertEquals(posStart + 3, feat.getStart());
        assertEquals(posStart + 7, feat.getEnd());

        feat = matchIter.next();
        assertEquals(posStart + 14, feat.getStart());
        assertEquals(posStart + 18, feat.getEnd());

        assertFalse(matchIter.hasNext());
    }

    @Test
    public void testConvertMotifToRegex_Basic() throws Exception{
        String motif = "ACTGACTGACTG";
        String regex = MotifFinderDialog.convertMotifToRegex(motif);
        assertEquals(motif, regex);
    }

    @Test
    public void testConvertMotifToRegex_02() throws Exception{
        String motif = "ACTGMACTGNACTSG";
        String regex = MotifFinderDialog.convertMotifToRegex(motif);

        assertTrue(regex.length() >= motif.length());
        assertEquals("ACTG[M,A,C]ACTG.ACT[S,G,C]G", regex);
    }

    @Test
    public void testExactSearch_EGFR() throws Exception{
        String motif = "CTTCGGGGAGCAGCGATGCGACCCTCCGGGACGGCCGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGCCCGGCGAGTCGGGCTCTGGAGGAAAAGAAAGGTAAGGGCGTGTCTCGCCGGCTCCCGCGCCGCCCCCGGATCGCGCCCCGGACCCCGCAGCCCGCCCAACCGCG";

        int expStart = 55054449;
        tstSearchGenome_EGFR(motif, expStart, expStart + motif.length());
    }

    @Test
    public void testAmbiguousSearch_EGFR() throws Exception{
        String motif = "CTTYKSVDAGCAGNGATGCRRCCCYCCGGGACGGCCGGGNCAGCGCKCCBGGCGCDGCTGGCTGCGCTCTGCCCGGCGAGTCGGGCTCTGGAGGRMWHGAAAGGNNVGGGCGTGTCTCGCCGGCTCCCGCGCCGCCCCCGGATCGCGCCCCGGACCCCGCAGCCCGCCCAACCGCG";

        int expStart = 55054449;
        String pattern = MotifFinderDialog.convertMotifToRegex(motif);
        tstSearchGenome_EGFR(pattern, expStart, expStart + motif.length());
    }

    @Test
    public void testSearchNoResult() throws Exception{
        tstSearchGenomeSingResult("GATCRYMKSWHBVDNGATCGATCGATCGATCGATCGATCGATCGATCGATC", "chr7",
                0, 100000, -1, -1);
    }

    /**
     * Test searching for a motif which has 2 overlapping hits
     * @throws Exception
     */
    @Test
    public void testSearchOverlapping() throws Exception{
        String motif = "ATGCATGCATGC";
        MotifFinderSource source = new MotifFinderSource(motif, genome);

        String queryChr = "chr8";
        int queryStart = 40823007;
        int queryEnd = 40863995;

        Iterator<Feature> iter = source.getFeatures(queryChr, queryStart, queryEnd);
        Feature feat = iter.next();

        int expFeatureStart = 40843500;
        int exFeatureEnd = 40843512;

        assertEquals(expFeatureStart, feat.getStart());
        assertEquals(exFeatureEnd, feat.getEnd());

        feat = iter.next();
        assertEquals(expFeatureStart + 4, feat.getStart());
        assertEquals(exFeatureEnd + 4, feat.getEnd());
    }

    /**
     * Test searching chromosome 1 for nonexistent motif. Test is timed at 30 seconds,
     * this is a loose performance test. Takes ~13 seconds on my machine
     * @throws Exception
     */
    @Test
    public void testSearchWholeChromo() throws Exception{
        String chromo = "chr1";
        String umotif = "GATCRYMKSWHBVDNGATCGATCGATCGATCGATCGATCGATCGATCGATC";
        int reps = 20;

        //Really long string, probably no matches
        String motif = umotif;
        for(int ii=0; ii < reps; ii++) motif += umotif;

        tstSearchGenomeSingResult(motif, chromo,
                0, genome.getChromosome(chromo).getLength(), -1, -1);
    }

    public void tstSearchGenome_EGFR(String motif, int expStart, int expEnd) throws Exception {
        //Our favorite, EGFR
        String chr = "chr7";
        int start = 55052219;
        int end = 55244525;

        tstSearchGenomeSingResult(motif, chr, start, end, expStart, expEnd);
    }

    /**
     * Test searching for pattern in the specified region, asserting that at most 1 feature
     * is found and starts at {@code expFeatureStart}. If expFeatureStart < 0, it is assumed
     * no features should be found
     * @param pattern
     * @param chr
     * @param start
     * @param end
     * @param expFeatStart
     * @throws Exception
     */
    public void tstSearchGenomeSingResult(String pattern, String chr, int start, int end, int expFeatStart, int expFeatureEnd) throws Exception {
        MotifFinderSource source = new MotifFinderSource(pattern, genome);

        Iterator<Feature> iter = source.getFeatures(chr, start, end);
        Feature feat = iter.next();

        if(expFeatStart >= 0){
            assertNotNull(feat);

            assertEquals(expFeatStart, feat.getStart());

            if(expFeatureEnd > 0){
                assertEquals(expFeatureEnd, feat.getEnd());
            }
        }else{
            assertNull(feat);
        }

        assertFalse(iter.hasNext());
    }
}
