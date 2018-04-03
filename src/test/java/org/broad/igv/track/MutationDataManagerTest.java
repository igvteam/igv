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

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.Mutation;
import org.broad.igv.util.ResourceLocator;
import org.junit.Test;

import java.util.Iterator;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 4/9/13
 *         Time: 10:01 AM
 */
public class MutationDataManagerTest extends AbstractHeadlessTest {

    String path = org.broad.igv.util.TestUtils.DATA_DIR + "maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf";

    @Test
    public void testGetFeatures() throws Exception {

        String sample = "TCGA-02-0055-01A-01W";
        String chr = "chr1";
        int start = 97575730;
        int end = 123482409;

        MutationFeatureSource.MutationDataManager mgr = new MutationFeatureSource.MutationDataManager(new ResourceLocator(path), genome);
        Iterator<Mutation> mutations =  mgr.getFeatures(sample, chr, start, end);

        int mutationCount = 0;
        while(mutations.hasNext()) {
            Mutation m = mutations.next();
            assertEquals(chr, m.getChr());

            if(m.getStart() >= start && m.getEnd() <= end) {
                mutationCount++;
            }
        }

        // There are 2 mutations for this sample in this interval
        assertEquals(2, mutationCount);
    }


}
