/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature.tribble;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 4/8/13
 *         Time: 2:51 PM
 */
public class MUTCodecTest extends AbstractHeadlessTest {

    String path = TestUtils.DATA_DIR + "maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf";



    @Test
    public void testIsMutationAnnotationFile() throws Exception {
        assertTrue(MUTCodec.isMutationAnnotationFile(path));
    }

    @Test
    public void testGetSamples() throws Exception {

        int expectedSampleCount = 146;
        MUTCodec codec = new MUTCodec(path, genome);
        String [] samples = codec.getSamples();
        assertEquals(expectedSampleCount, samples.length);
        assertTrue(samples[0].startsWith("TCGA"));

    }


}
