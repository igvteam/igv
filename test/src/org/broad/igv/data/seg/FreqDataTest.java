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

package org.broad.igv.data.seg;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;

/**
 * @author jrobinso
 * @date Oct 13, 2010
 */
public class FreqDataTest extends AbstractHeadlessTest {

    @Test
    public void test() throws IOException {

        Genome genome = TestUtils.loadGenome();

        String segfile = TestUtils.DATA_DIR + "seg/canFam2_hg18.seg";


        SegmentedDataSet sd = new SegmentedAsciiDataSet(new ResourceLocator(segfile), genome);

        FreqData fd = new FreqData(sd, genome);
    }
}
