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

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.*;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class AlignmentPackerTest extends AbstractHeadlessTest {

    String path = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
    String chr = "chr1";
    int start = 151666494;
    int end = start + 1000;
    boolean contained = false;


    private AlignmentInterval getAlignmentInterval() throws Exception {
        ResourceLocator rl = new ResourceLocator(path);
        AlignmentReader samReader = AlignmentReaderFactory.getReader(rl);
        CloseableIterator<Alignment> iter = samReader.query(chr, start, end, contained);
        List<Alignment> list = new ArrayList<Alignment>();
        while(iter.hasNext()){
            list.add(iter.next());
        }
        AlignmentInterval interval = new AlignmentInterval(chr, start, end, list, null, null, null);
        return interval;
    }

    /**
     * Test of packAlignments method, of class AlignmentPacker.
     */
    @Test
    public void testPackAlignments() throws Exception {

        ///////////////////////////
        /*
        boolean showDuplicates = false;
        int qualityThreshold = 0;
        int maxLevels = 1000;
        */
        AlignmentInterval interval = getAlignmentInterval();

        Map<String, List<Row>> result = (new AlignmentPacker()).packAlignments(interval, new AlignmentTrack.RenderOptions());
        assertEquals(1, result.size());
        for (List<Row> alignmentrows : result.values()) {
            for (Row alignmentrow : alignmentrows) {
                List<Alignment> alignments = alignmentrow.alignments;
                for (int ii = 1; ii < alignments.size(); ii++) {
                    assertTrue(alignments.get(ii).getAlignmentStart() > alignments.get(ii - 1).getAlignmentStart());
                    assertTrue(alignments.get(ii).getAlignmentStart() - alignments.get(ii - 1).getAlignmentEnd() >= AlignmentPacker.MIN_ALIGNMENT_SPACING);
                }
            }
        }


    }

    @Test
    public void testGroupAlignmentsPairOrientation() throws Exception {
        int expSize = 3; //AlignmentTrack.OrientationType.values().length;
        Map<String, List<Row>> result = tstGroupAlignments(AlignmentTrack.GroupOption.PAIR_ORIENTATION, expSize);
    }

    public Map<String, List<Row>> tstGroupAlignments(AlignmentTrack.GroupOption groupOption, int expSize) throws Exception {

        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        renderOptions.groupByOption = groupOption;

        AlignmentInterval interval = getAlignmentInterval();
        Map<String, List<Row>> result = (new AlignmentPacker()).packAlignments(interval, renderOptions);
        Set<String> names = result.keySet();
        //names.removeAll(Arrays.asList("", null));

        assertEquals(expSize, names.size());
        return result;

    }


}