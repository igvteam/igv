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

package org.broad.igv.methyl;

import org.broad.igv.feature.genome.Genome;
import org.junit.Test;

import java.util.Iterator;

import static junit.framework.Assert.assertTrue;

/**
 * Test of data source for Michael Ziller's hacked up bigbed format for methylation scores.
 *
 * @author Jim Robinson
 * @date 4/19/12
 */
public class BBMethylDataSourceTest {

    @Test
    public void pointlessTest() {
        //Must have at least 1 test method
    }

    //@Test
    public void testQuery() throws Exception {

        String testFile = "path-to-test-file.bb";
        String chr = "chr7";
        int start = 26133475;
        int end = 27333475;
        Genome genome = null;

        MethylDataSource source = new BBMethylDataSource(testFile, genome);
        Iterator<MethylScore> iter = source.query(chr, start, end);

        // Should find at least 1 score.
        assertTrue(iter.hasNext());

        while (iter.hasNext()) {
            MethylScore score = iter.next();
            assertTrue(score.getEnd() >= start && score.getStart() <= end);
            assertTrue(score.getScore() >= 0 && score.getScore() <= 100);
        }
    }
}
