package org.igv.ui.panel;

import org.igv.lists.GeneList;
import org.igv.ui.AbstractHeadedTest;
import org.igv.ui.IGV;
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
        GeneList geneList = new GeneList("", loci);

        IGV.getInstance().getSession().setCurrentGeneList(geneList);
        IGV.getInstance().resetFrames();
        IGV.getInstance().waitForNotify(5000);

        ReferenceFrameTest.assertLociFramesConsistent(loci, FrameManager.getFrames());
    }

}
