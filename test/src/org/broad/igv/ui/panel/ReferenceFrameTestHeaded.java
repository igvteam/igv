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
