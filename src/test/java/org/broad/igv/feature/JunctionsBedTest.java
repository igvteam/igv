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
