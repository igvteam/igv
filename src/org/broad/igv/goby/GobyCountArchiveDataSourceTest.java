package org.broad.igv.goby;

import org.broad.igv.feature.LocusScore;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * Test that the count goby data source can load data.
 * @author Fabien Campagne
 *         Date: 6/12/11
 *         Time: 1:34 PM
 */
public class GobyCountArchiveDataSourceTest {
    @Test
    public void iteratePickrell() {

        GobyCountArchiveDataSource counts = new GobyCountArchiveDataSource(new File("test/data/goby/counts/CIIIBUD-Pickrell-2010NA19204_argonne.counts"));
        assertTrue(counts.getAvailableWindowFunctions().size() > 0);
        List<LocusScore> list1 = counts.getSummaryScoresForRange("chr1", 10000, 2000000, 1);
        assertTrue(list1.size() > 0);
        List<LocusScore> list2 = counts.getSummaryScoresForRange("1", 10000, 2000000, 1);
        assertTrue(list2.size() > 0);
        assertEquals(list1.size(), list2.size());
    }
}
