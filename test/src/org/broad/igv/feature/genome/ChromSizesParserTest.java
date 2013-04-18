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

package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 4/16/13
 *         Time: 2:43 PM
 */
public class ChromSizesParserTest {

    String testFile = TestUtils.DATA_DIR + "genomes/hg19.chrom.sizes";

    @Test
    public void testParse() throws Exception {

        int expectedCount = 24;

        List<Chromosome> chromosomes = ChromSizesParser.parse(testFile);
        assertEquals(expectedCount, chromosomes.size());

        // chr10	135534747	/gbdb/hg19/hg19.2bit
        Chromosome chr10 = chromosomes.get(9);
        assertEquals("chr10", chr10.getName());
        assertEquals(135534747, chr10.getLength());
    }
}
