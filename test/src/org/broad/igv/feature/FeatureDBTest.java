/*
 * Copyright (c) 2007-2012 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.feature;

import org.broad.igv.Globals;
import org.broad.igv.TestInformation;
import org.broad.igv.tools.IgvTools;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.Map;

import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2011/12/15
 */
public class FeatureDBTest {

    static String dataFileName = TestInformation.DATA_DIR + "/genomes/hg18.genome";

    @Before
    public void setUp() throws Exception {
        Globals.setHeadless(true);
        IgvTools.loadGenome(dataFileName, true);
    }

    @After
    public void tearDown() {
        FeatureDB.clearFeatures();
    }

    @Test
    public void testPartialMatch() throws Exception {
        String chkstr = "ABC";
        Map<String, NamedFeature> fMap = FeatureDB.getFeatures(chkstr);

        for (String k : fMap.keySet()) {

            assertTrue(k.startsWith(chkstr));
        }

    }


}
