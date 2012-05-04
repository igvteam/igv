/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature;

import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
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
