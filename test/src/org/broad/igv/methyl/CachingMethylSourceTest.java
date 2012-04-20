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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 4/19/12
 */
public class CachingMethylSourceTest {
   // @Test
    public void testQuery() throws Exception {

        String testFile = "path-to-test-file.bb";
        String chr = "chr7";
        int start = 26133475;
        int end = 27333475;
        Genome genome = null;

        MethylDataSource source = new BBMethylDataSource(testFile, genome);
        Iterator<MethylScore> iter = source.query(chr, start, end);
        List<MethylScore> nonCachedScores = new ArrayList<MethylScore>();
        while (iter.hasNext()) {
            nonCachedScores.add(iter.next());
        }

        int maxTileCount = 100000;
        int tileSize = 100;
        MethylDataSource cachedSource = new CachingMethylSource( new BBMethylDataSource(testFile, genome),
                maxTileCount, tileSize);
        Iterator<MethylScore> iter2 = cachedSource.query(chr, start, end);
        List<MethylScore> cachedScores = new ArrayList<MethylScore>();
        while (iter2.hasNext()) {
            cachedScores.add(iter2.next());
        }

        assertTrue(cachedScores.size() > 0);
        assertEquals(cachedScores.size(), nonCachedScores.size());
        for(int i=0; i<cachedScores.size(); i++) {
            MethylScore cScore = cachedScores.get(i);
            MethylScore score = nonCachedScores.get(i);
            assertEquals(cScore.getStart(), score.getStart());
            assertEquals(cScore.getChr(), score.getChr());
        }
    }
}
