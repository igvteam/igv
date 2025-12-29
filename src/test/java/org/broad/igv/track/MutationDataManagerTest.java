package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.Mutation;
import org.broad.igv.util.ResourceLocator;
import org.junit.Test;

import java.util.Iterator;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 4/9/13
 *         Time: 10:01 AM
 */
public class MutationDataManagerTest extends AbstractHeadlessTest {

    String path = org.broad.igv.util.TestUtils.DATA_DIR + "maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf";

    @Test
    public void testGetFeatures() throws Exception {

        String sample = "TCGA-02-0055-01A-01W";
        String chr = "chr1";
        int start = 97575730;
        int end = 123482409;

        MutationFeatureSource.MutationDataManager mgr = new MutationFeatureSource.MutationDataManager(new ResourceLocator(path), genome);
        Iterator<Mutation> mutations =  mgr.getFeatures(sample, chr, start, end);

        int mutationCount = 0;
        while(mutations.hasNext()) {
            Mutation m = mutations.next();
            assertEquals(chr, m.getChr());

            if(m.getStart() >= start && m.getEnd() <= end) {
                mutationCount++;
            }
        }

        // There are 2 mutations for this sample in this interval
        assertEquals(2, mutationCount);
    }


}
