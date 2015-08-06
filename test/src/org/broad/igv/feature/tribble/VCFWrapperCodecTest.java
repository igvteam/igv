/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.feature.tribble;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.vcf.VCFVariant;
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
