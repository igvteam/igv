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

package org.broad.igv.tools.motiffinder;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Strand;
import htsjdk.tribble.Feature;
import org.junit.Test;

import java.util.Iterator;
import java.util.regex.Matcher;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2013-Jan-22
 */
public class MotifFinderSourceTest extends AbstractHeadlessTest{


    String shortSeq = "GACTCTGACTGGAGTCTGATCAG";
    //revComp =       "CTGATCAGACTCCAGTCAGAGTC";

    @Test
    public void testSearchShortSeqPosStrand() throws Exception{
        tstSearchShortSeq(Strand.POSITIVE, 2);
    }

    @Test
    public void testSearchShortSeqNegStrand() throws Exception{
        tstSearchShortSeq(Strand.NEGATIVE, 1);
    }

    private void tstSearchShortSeq(Strand strand, int expFeatures) throws Exception{
        int posStart = 0;//19389;
        //Motif is twice on positive strand, once on negative
        String motif = "ACT";

        Iterator<Feature> matchIter = MotifFinderSource.search(motif, strand, null, posStart, shortSeq.getBytes());
        assertTrue("No matches found for motif " + motif, matchIter.hasNext());

        int count = 0;
        int lastStart = -1;
        while(matchIter.hasNext()){
            Feature feat = matchIter.next();
            count++;
            assertTrue("Features out of order", feat.getStart() >= lastStart);
            checkPatternMatches(motif, feat, shortSeq.getBytes());
            lastStart = feat.getStart();
        }
        assertEquals("Unexpected number of features", expFeatures, count);
    }

    @Test
    public void testBasicSearchNegStrand_end() throws Exception{
        int posStart = 19389;
        //Motif is at the end of this sequence on the negative strand
        String motif = "CTGATC";
        Strand strand = Strand.NEGATIVE;

        Iterator<Feature> matchIter = MotifFinderSource.searchSingleStrand(motif, strand, null, posStart, shortSeq.getBytes());
        assertTrue("No matches found for motif " + motif, matchIter.hasNext());
        Feature feat = matchIter.next();

        int expStart = posStart + shortSeq.length() - motif.length();
        int expEnd = expStart + motif.length();
        assertEquals(expStart, feat.getStart());
        assertEquals(expEnd, feat.getEnd());
        checkPatternMatches(motif, feat, shortSeq.getBytes());

        assertFalse(matchIter.hasNext());
    }

    @Test
    public void testBasicSearchNegStrand_few() throws Exception{
        int posStart = 19389;
        //Should have several hits
        String motif = "TG";
        String revSeq = "CA";
        Strand strand = Strand.NEGATIVE;

        Iterator<Feature> matchIter = MotifFinderSource.searchSingleStrand(motif, strand, null, posStart, shortSeq.getBytes());
        assertTrue("No matches found for motif " + motif, matchIter.hasNext());

        while(matchIter.hasNext()){
            Feature feat = matchIter.next();
            checkPatternMatches(motif, feat, shortSeq.getBytes());
        }
    }

    @Test
    public void testSearchSingleStrandPos() throws Exception{
        int posStart = 19389;
        String motif = "TCTG";
        Strand strand = Strand.POSITIVE;

        Iterator<Feature> matchIter = MotifFinderSource.searchSingleStrand(motif, strand, null, posStart, shortSeq.getBytes());
        assertTrue("No matches found for motif " + motif, matchIter.hasNext());
        Feature feat = matchIter.next();

        assertEquals(posStart + 3, feat.getStart());
        assertEquals(posStart + 7, feat.getEnd());
        checkPatternMatches(motif, feat, shortSeq.getBytes());

        feat = matchIter.next();
        assertEquals(posStart + 14, feat.getStart());
        assertEquals(posStart + 18, feat.getEnd());
        checkPatternMatches(motif, feat, shortSeq.getBytes());


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
        tstSearchGenomeSingResult("GATCRYMKSWHBVDNGATCGATCGATCGATCGATCGATCGATCGATCGATC", Strand.POSITIVE, "chr7",
                0, 100000, -1, -1);
    }

    /**
     * Test searching for a motif which has 2 overlapping hits
     * @throws Exception
     */
    @Test
    public void testSearchOverlapping() throws Exception{
        String motif = "ATGCATGCATGC";
        MotifFinderSource source = new MotifFinderSource(motif, Strand.POSITIVE, genome);

        String queryChr = "chr8";
        int queryStart = 40823007;
        int queryEnd = 40863995;

        Iterator<Feature> iter = source.getFeatures(queryChr, queryStart, queryEnd);
        BasicFeature feat = (BasicFeature) iter.next();

        int expFeatureStart = 40843500;
        int exFeatureEnd = 40843512;

        assertEquals(expFeatureStart, feat.getStart());
        assertEquals(exFeatureEnd, feat.getEnd());
        checkPatternMatches(motif, feat);

        //Negative strand hit
//        feat = (BasicFeature) iter.next();
//        assertEquals(expFeatureStart + 2, feat.getStart());
//        assertEquals(exFeatureEnd + 2, feat.getEnd());
//        checkPatternMatches(motif, feat);

        feat = (BasicFeature) iter.next();
        assertEquals(expFeatureStart + 4, feat.getStart());
        assertEquals(exFeatureEnd + 4, feat.getEnd());
        checkPatternMatches(motif, feat);
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

        tstSearchGenomeSingResult(motif, Strand.POSITIVE, chromo,
                0, genome.getChromosome(chromo).getLength(), -1, -1);
    }

    private void tstSearchGenome_EGFR(String motif, int expStart, int expEnd) throws Exception {
        //Our favorite, EGFR
        String chr = "chr7";
        int start = 55052219;
        int end = 55244525;

        tstSearchGenomeSingResult(motif, Strand.POSITIVE, chr, start, end, expStart, expEnd);
    }

    /**
     * Test searching for pattern in the specified region, asserting that at most 1 feature
     * is found and starts at {@code expFeatureStart}. If expFeatureStart < 0, it is assumed
     * no features should be found
     * @param pattern
     * @param strand
     * @param chr
     * @param start
     * @param end
     * @param expFeatStart
     * @throws Exception
     */
    private void tstSearchGenomeSingResult(String pattern, Strand strand, String chr, int start, int end, int expFeatStart, int expFeatureEnd) throws Exception {
        MotifFinderSource source = new MotifFinderSource(pattern, strand, genome);

        Iterator<Feature> iter = source.getFeatures(chr, start, end);

        if(expFeatStart >= 0){

            Feature feat = iter.next();
            assertNotNull(feat);

            assertEquals(expFeatStart, feat.getStart());

            if(expFeatureEnd > 0){
                assertEquals(expFeatureEnd, feat.getEnd());
            }

            checkPatternMatches(pattern, feat);
        }else{
            assertFalse("Iterator should be empty but isn't", iter.hasNext());
        }

        assertFalse(iter.hasNext());
    }

    private void checkPatternMatches(String pattern, Feature feature){
        byte[] sequence = genome.getSequence(feature.getChr(), feature.getStart(), feature.getEnd());
        checkPatternMatches(pattern, feature, sequence);
    }

    private void checkPatternMatches(String pattern, Feature feature, byte[] sequence){
        Matcher matcher = MotifFinderSource.getMatcher(pattern, ((BasicFeature) feature).getStrand(), sequence);
        assertTrue("No match found for pattern " + pattern + " at " + feature.getStart(), matcher.find());
    }
}
