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
import org.broad.igv.util.TestUtils;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.junit.Test;

import java.io.File;
import java.util.Iterator;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 * @date Aug 9, 2010
 */
public class TribbleIndexTest {

    IgvTools igvTools = new IgvTools();

    @Test
    /**
     * Compare linear and interval tree index
     */
    public void testIntervalTree() throws Exception {
        //chr2:179,222,066-179,262,059<- CONTAINS TTN
        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";
        String chr = "chr2";
        int start = 179222066;
        int end = 179262059;
        int expectedCount = 3;

        // Linear index
        File indexFile = new File(bedFile + ".idx");
        if (indexFile.exists()) {
            indexFile.delete();
        }
        igvTools.doIndex(bedFile, IgvTools.LINEAR_INDEX, 1000);
        indexFile.deleteOnExit();


        BasicFeatureSource bfr = BasicFeatureSource.getFeatureSource(bedFile, new IGVBEDCodec());
        Iterator<BasicFeature> iter = bfr.query(chr, start, end);
        int count = 0;
        while (iter.hasNext()) {
            BasicFeature feat = iter.next();
            count++;
        }


        // Interval index
        if (indexFile.exists()) {
            indexFile.delete();
        }
        igvTools.doIndex(bedFile, IgvTools.INTERVAL_INDEX, 5);
        indexFile.deleteOnExit();


        bfr = BasicFeatureSource.getFeatureSource(bedFile, new IGVBEDCodec());
        iter = bfr.query(chr, start, end);
        int countInterval = 0;
        while (iter.hasNext()) {
            BasicFeature feat = iter.next();
            String tmp = "" + feat.getStart() + " " + feat.getEnd();
            assertTrue(tmp, feat.getStart() <= end && feat.getEnd() >= (start - 1));
            countInterval++;
        }

        assertEquals(count, countInterval);
        assertEquals(expectedCount, count);

    }

    @Test
    public void testReadSingleVCF() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "/vcf/indel_variants_onerow.vcf";
        String chr = "chr9";
        // Linear index
        File indexFile = new File(bedFile + ".idx");
        if (indexFile.exists()) {
            indexFile.delete();
        }
        igvTools.doIndex(bedFile, IgvTools.LINEAR_INDEX, 1000);
        indexFile.deleteOnExit();


        BasicFeatureSource bfr = BasicFeatureSource.getFeatureSource(bedFile, new VCFCodec());
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

    }

}
