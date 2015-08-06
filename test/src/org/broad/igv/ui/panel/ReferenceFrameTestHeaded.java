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

package org.broad.igv.ui.panel;

import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

/**
 * User: jacob
 * Date: 2013-Mar-20
 */
public class ReferenceFrameTestHeaded extends AbstractHeadedTest{


    @Override
    public void setUp() throws Exception {
        super.setUp();
        ReferenceFrameTest.RFTSetup();
    }

    /**
     * Test that when we view a gene list it displays correctly
     */
    @Test
    public void testGeneListView(){
        List<String> loci = Arrays.asList("chr1:10000-20000", "chr2:30000-40000");
        GeneList geneList = new GeneList("", loci, false);

        IGV.getInstance().getSession().setCurrentGeneList(geneList);
        IGV.getInstance().resetFrames();
        IGV.getInstance().waitForNotify(5000);

        ReferenceFrameTest.assertLociFramesConsistent(loci, FrameManager.getFrames());
    }

}
