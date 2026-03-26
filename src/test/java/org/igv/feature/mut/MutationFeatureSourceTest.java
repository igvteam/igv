package org.igv.feature.mut;

import org.igv.AbstractHeadlessTest;
import org.igv.feature.LocusScore;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Oct-23
 */
public class MutationFeatureSourceTest extends AbstractHeadlessTest {
    @Test

    public void testMAF() throws Exception {

        String testPath = TestUtils.DATA_DIR + "maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf.gz";

        MutationFeatureSource track  = new MutationFeatureSource(new ResourceLocator(testPath), genome);

        String testSampleId = "TCGA-02-0024-01B-01W";
        String testChr = "chr12";
        int testStart = 9150104 - 1;
        int testEnd = testStart + 1;

        //Look for the feature at line 4 of the file:
        //A2M	2	genome.wustl.edu	36	12	9150104	9150104	+	Missense_Mutation	SNP	C	C	G	Unknown	Unknown	TCGA-02-0024-01B-01W	TCGA-02-0024-10A-01W	C	C	C	G	C	C	Valid	Valid	Somatic	Phase_I

        List<LocusScore> features = track.getSegments(testSampleId, testChr);

        boolean found = false;
       for(LocusScore feature : features) {
            Mutation mut = (Mutation) feature;
            if (mut.getSampleId().equals(testSampleId)) {
                assertEquals(testChr, mut.getChr());
                if(testStart == mut.getStart()) {
                    assertEquals(testStart, mut.getStart());
                    assertEquals(testEnd, mut.getEnd());
                    assertEquals("Missense_Mutation", mut.getMutationType());
                    assertEquals("CCG", mut.refAllele + mut.altAllele1 + mut.altAllele2);
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            org.junit.Assert.fail("Expected mutation not found in features");
        }

        List<String> samples = track.getSampleNames();

        assertEquals(145, samples.size());

    }
}
