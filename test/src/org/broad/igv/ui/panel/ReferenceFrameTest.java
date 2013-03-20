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

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Locus;
import org.broad.igv.lists.GeneList;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2013-Mar-20
 */
public class ReferenceFrameTest extends AbstractHeadlessTest{

    private ReferenceFrame frame;

    @Override
    public void setUp() throws Exception {
        super.setUp();
        RFTSetup();
        frame = new ReferenceFrame("testFrame");
        frame.setChromosomeName(Globals.CHR_ALL);
        frame.setBounds(0, 500);
    }

    @Override
    public void tearDown() throws Exception {
        super.tearDown();
        frame = null;
    }

    @Test
    public void testSetBounds_chr1() throws Exception{
        frame.setChromosomeName("chr1");

        double oldLocScale = frame.locationScale;
        double oldOrigin = frame.getOrigin();
        double oldEnd = frame.getEnd();
        int oldWidth = frame.getWidthInPixels();

        double multiplier = 1.2d;
        int newWidth = (int) (multiplier * oldWidth);
        frame.setBounds(0, newWidth);

        assertEquals(oldOrigin, frame.getOrigin());
        assertEquals(oldEnd, frame.getEnd(), 1.0);
        assertEquals(oldLocScale, frame.locationScale * multiplier, 2.0);

        assertConsistent();
    }

    @Test
    public void testJumpTo_00() throws Exception{
        Locus locus = new Locus("chr1", 1000, 2000);
        frame.jumpTo(locus);

        assertEquals(locus.getChr(), frame.getChrName());
        assertEquals(locus.getStart(), frame.getOrigin(), 1.0);
        assertEquals(locus.getEnd(), frame.getEnd(), 1.0);
        assertConsistent();
    }

    @Test
    public void testJumpTo_01() throws Exception{
        Locus locus = new Locus("chr1", 1000, 2000);
        frame.jumpTo(locus);
        double oldLocScale = frame.getScale();
        int oldZoom = frame.getZoom();
        int oldMidPoint = frame.getMidpoint();
        double oldCenter = frame.getCenter();

        int delta = 12344;
        frame.jumpTo(locus.getChr(), locus.getStart() + delta, locus.getEnd() + delta);

        assertEquals(oldLocScale, frame.getScale());
        assertEquals(oldZoom, frame.getZoom());
        assertEquals(oldMidPoint, frame.getMidpoint());

        assertEquals(oldCenter + delta, frame.getCenter(), 1.0);

        assertConsistent();
    }

    @Test
    public void testDoZoomIncrement() throws Exception{
        frame.setChromosomeName("chr1");

        double oldCenter = frame.getCenter();
        double oldLocScale = frame.getScale();
        int oldZoom = frame.getZoom();
        int incrAmount = 2;

        frame.doZoomIncrement(incrAmount);

        assertConsistent();

        assertEquals(oldZoom + incrAmount, frame.zoom);
        assertEquals(oldCenter, frame.getCenter(), 2.0);

        double scaleMultiplier = Math.pow(2.0, incrAmount);
        assertEquals(oldLocScale, scaleMultiplier * frame.getScale(), 1.0);
    }

    /**
     * Test that when we view a gene list it displays correctly
     */
    @Test
    public void testGeneListView(){

        int paneSize = 300;

        List<String> loci = Arrays.asList("chr1:10000-20000", "chr2:30000-40000");
        GeneList geneList = new GeneList("", loci, false);
        FrameManager.resetFrames(geneList);
        List<ReferenceFrame> frameList = FrameManager.getFrames();

        //Simulate doing layout
        for(int ff = 0; ff < loci.size(); ff++){
            frameList.get(ff).setBounds(ff*paneSize, paneSize);
        }

        assertLociFramesConsistent(loci, frameList);
    }

    private void assertConsistent(){
        assertConsistent(frame);
    }

    static void assertConsistent(ReferenceFrame frame){
        assertEquals(frame.getEnd(), frame.origin + frame.locationScale * frame.widthInPixels);
        assertEquals(frame.getnTiles(), Math.pow(2, frame.zoom));
    }

    static void assertLociFramesConsistent(List<String> loci, List<ReferenceFrame> frameList){
        assertEquals(loci.size(), frameList.size());

        for(int ii=0; ii < loci.size(); ii++){
            ReferenceFrame frame = frameList.get(ii);
            ReferenceFrameTest.assertConsistent(frame);

            Locus locus = Locus.fromString(loci.get(ii));

            assertEquals(locus.getChr(), frame.getChrName());
            assertEquals(locus.getStart() - 1, frame.getOrigin(), 0.5);
            assertEquals(locus.getEnd(), frame.getEnd(), 0.5);
        }
    }

    public static void RFTSetup() {
        PreferenceManager.getInstance().put(PreferenceManager.FLANKING_REGION, "" + 0);
    }
}
