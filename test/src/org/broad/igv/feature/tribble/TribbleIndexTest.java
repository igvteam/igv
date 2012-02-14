/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.feature.tribble;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.junit.Assert;
import org.junit.Test;

import java.util.*;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 * @date Aug 9, 2010
 */
public class TribbleIndexTest {

    IgvTools igvTools = new IgvTools();


    /**
     * chr2	1	200000000	LONG_FEATURE
     * ...
     * chr2	179098961	179380395	Hs.134602
     * chr2	179209546	179287210	Hs.620337
     * chr2	179266309	179266748	Hs.609465
     * chr2	179296428	179300012	Hs.623987
     * chr2	179302952	179303488	Hs.594545
     */

    @Test
    public void testLinearIndex() throws Exception {

        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";
        String chr = "chr2";
        int start = 179266309 - 1;
        int end = 179303488 + 1;
        int expectedCount = 6;

        Set<String> expectedNames = new HashSet<String>(Arrays.asList("Hs.134602", "Hs.620337", "Hs.609465", "Hs.623987",
                "Hs.594545", "LONG_FEATURE"));

        // Interval index
        TestUtils.createIndex(bedFile, IgvTools.LINEAR_INDEX, 500);

        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, new IGVBEDCodec());
        Iterator<BasicFeature> iter = bfr.query(chr, start, end);
        int countInterval = 0;
        while (iter.hasNext()) {
            BasicFeature feature = iter.next();
            Assert.assertTrue(feature.getEnd() >= start && feature.getStart() <= end);
            Assert.assertTrue(expectedNames.contains(feature.getName()));
            countInterval++;
        }

        assertEquals(expectedCount, countInterval);

    }

    @Test
    /**
     * Test interval tree index
     * chr2	179098961	179380395	Hs.134602
     * chr2	179209546	179287210	Hs.620337
     * chr2	179266309	179266748	Hs.609465
     * chr2	179296428	179300012	Hs.623987
     * chr2	179302952	179303488	Hs.594545
     *
     */
    public void testIntervalTree() throws Exception {
        //chr2:179,222,066-179,262,059<- CONTAINS TTN
        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";
        String chr = "chr2";
        int start = 179266309 - 1;
        int end = 179303488 + 1;
        int expectedCount = 6;

        Set<String> expectedNames = new HashSet<String>(Arrays.asList("Hs.134602", "Hs.620337", "Hs.609465", "Hs.623987",
                "Hs.594545", "LONG_FEATURE"));

        // Interval index
        TestUtils.createIndex(bedFile, IgvTools.INTERVAL_INDEX, 1);

        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, new IGVBEDCodec());
        Iterator<BasicFeature> iter = bfr.query(chr, start, end);
        int countInterval = 0;
        while (iter.hasNext()) {
            BasicFeature feature = iter.next();
            Assert.assertTrue(feature.getEnd() >= start && feature.getStart() <= end);
            Assert.assertTrue(expectedNames.contains(feature.getName()));
            countInterval++;
        }

        assertEquals(expectedCount, countInterval);
    }

    @Test
    public void testReadSingleVCF() throws Exception {
        String file = TestUtils.DATA_DIR + "/vcf/indel_variants_onerow.vcf";
        String chr = "chr9";
        // Linear index
        TestUtils.createIndex(file);

        // First test query
        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(file, new VCFCodec());
        Iterator<org.broadinstitute.sting.utils.variantcontext.VariantContext> iter = bfr.query(chr, 5073767 - 5, 5073767 + 5);
        int count = 0;
        while (iter.hasNext()) {
            org.broadinstitute.sting.utils.variantcontext.VariantContext feat = iter.next();
            assertEquals("chr9", feat.getChr());
            assertEquals(feat.getStart(), 5073767);
            assertTrue(feat.hasAttribute("MapQs"));
            count++;
        }
        assertEquals(1, count);

        // Test non-indexed access (iterator)
        iter = bfr.iterator();
        count = 0;
        while (iter.hasNext()) {
            org.broadinstitute.sting.utils.variantcontext.VariantContext feat = iter.next();
            assertEquals("chr9", feat.getChr());
            assertEquals(feat.getStart(), 5073767);
            assertTrue(feat.hasAttribute("MapQs"));
            count++;
        }
        assertEquals(1, count);

        //Do similar as above, but have a different test file
        file = TestUtils.DATA_DIR + "/vcf/outputPileup.flt1.vcf";
        chr = "1";
        // Linear index
        TestUtils.createIndex(file);

        bfr = AbstractFeatureReader.getFeatureReader(file, new VCFCodec());
        iter = bfr.query(chr, 984163 - 5, 984163 + 5);
        count = 0;
        while (iter.hasNext()) {
            org.broadinstitute.sting.utils.variantcontext.VariantContext feat = iter.next();
            assertEquals(chr, feat.getChr());
            if (count == 0) {
                assertEquals(984163, feat.getStart());
            }
            count++;
        }
        assertEquals(1, count);

    }

    //@Test
    public void testReadVCFGui() throws Exception {
        IGV igv = TestUtils.startGUI();

        String file = TestUtils.DATA_DIR + "/vcf/outputPileup.flt1.vcf";
        TestUtils.createIndex(file);
        ResourceLocator locator = new ResourceLocator(file);
        //For files with 1 record, this threw a null pointer exception prior to r1595
        List<Track> tracks = igv.load(locator);

        assertEquals(1, tracks.size());

        //todo FIND a way to do this. We don't actually need/want to close IGV properly, just kill it
        //Also this doesn't work
        IGV.getMainFrame().setVisible(false);
        IGV.getMainFrame().dispose();

        //IGV.getRootPane().setVisible(false);
        //IGV.getRootPane().dispatchEvent(new WindowEvent(IGV.getMainFrame(), WindowEvent.WINDOW_CLOSING));
        //igv.doExitApplication();
        //System.out.println(IGV.getMainFrame());
    }

}
