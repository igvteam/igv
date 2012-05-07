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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.feature;

import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertNotNull;

/**
 * @author jrobinso
 */
public class TestLoadGeneManager {

    @Test
    public void main() throws IOException {
        TestUtils.setUpHeadless();
        IgvTools.loadGenome(TestUtils.DATA_DIR + "genomes/hg18.unittest.genome", true);
        Feature feature = FeatureDB.getFeature("EGFR");
        assertNotNull(feature);
    }

}
