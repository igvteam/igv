package org.igv.feature.tribble;

import org.igv.AbstractHeadlessTest;
import org.igv.track.TribbleFeatureSource;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.igv.variant.vcf.VCFVariant;
import org.junit.Test;

import java.util.Iterator;

import static org.junit.Assert.*;

/**
 * @author jacob
 * @date 2013-Jun-21
 */
public class VCFWrapperCodecTest extends AbstractHeadlessTest {

    /**
     * It is apparently a matter of some contention whether having a missing
     * field within a comma-separated list of fields should be legal VCF.
     * We've decided that IGV will accept these files, but since we use the picard VCF codec (and they don't want to)
     * we need to work around it.
     * @throws Exception
     */
    @Test
    public void testMissingFieldInCommaSeparated() throws Exception{

        String filePath = TestUtils.DATA_DIR + "vcf/missingFields.vcf";
        TestUtils.createIndex(filePath);

        TribbleFeatureSource src = TribbleFeatureSource.getFeatureSource(new ResourceLocator(filePath), genome);

        Iterator iter = src.getFeatures("chr2", 3321000, 13346000);

        int count = 0;
        boolean found = false;

        while(iter.hasNext()){
            VCFVariant vcfVariant = (VCFVariant) iter.next();
            assertNotNull(vcfVariant);

            if(vcfVariant.getStart() == 3796932){
                assertEquals(9, vcfVariant.getSampleNames().size());
                found = true;
            }

            count++;
        }

        assertEquals(26, count);
        assertTrue("Feature of interest not found", found);


    }
}
