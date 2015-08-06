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

package org.broad.igv.feature;

import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import org.junit.AfterClass;

import static org.junit.Assert.*;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Iterator;

/**
 * @author Jim Robinson
 * @date 4/26/12
 */
public class JunctionsBedTest {


    /**
     * Purpose of this test is to insure that a cufflinks "junctions.bed" file is parsed as a junction file, as
     * opposed to a plain bed file.
     *
     * @throws Exception
     */
    @Test
    public void testJunctionFile() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "bed/mini.junctions.bed";
        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, new IGVBEDCodec(), false);
        Iterator<BasicFeature> iter = bfr.iterator();
        while (iter.hasNext()) {
            BasicFeature feature = iter.next();
            assertTrue( feature instanceof SpliceJunctionFeature);
            SpliceJunctionFeature sjf = (SpliceJunctionFeature) feature;
            assertTrue(sjf.getJunctionDepth() > 0);
        }
    }


}
