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

package org.broad.igv.variant.vcf;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * The testXXXStart methods exist because we offset the start position based on the
 * bases common between reference and alternate alleles.
 * @see org.broad.igv.variant.vcf.VCFVariant#calcStart()
 *
 *
 * @author jacob
 * @date 2013-Apr-03
 */
public class VCFVariantTest extends AbstractHeadlessTest {

    private VCFVariant get_hc_mod_Feat(String chr, int start, int end) throws Exception{
        String filePath = TestUtils.DATA_DIR + "vcf/hc_mod.vcf";
        TestUtils.createIndex(filePath);
        TribbleFeatureSource src = TribbleFeatureSource.getFeatureSource(new ResourceLocator(filePath), genome);

        return (VCFVariant) (src.getFeatures(chr, start, end)).next();
    }

    @Test
    public void testSNPStart() throws Exception{
        int fileFeatStart = 10499005;
        tstStart("chr20", fileFeatStart, VariantContext.Type.SNP, fileFeatStart - 1);
    }

    @Test
    public void testBlockSubStart() throws Exception{
        int fileFeatStart = 10499073;
        tstStart("chr20", fileFeatStart, VariantContext.Type.MNP, fileFeatStart - 1);
    }

    @Test
    public void testDELStart_pref() throws Exception{
        int fileFeatStart = 10499099;
        tstStart("chr20", fileFeatStart, VariantContext.Type.INDEL, fileFeatStart - 1 + 1);
    }

    @Test
    public void testDELStart_nopref() throws Exception{
        int fileFeatStart = 10499291;
        tstStart("chr20", fileFeatStart, VariantContext.Type.INDEL, fileFeatStart - 1);
    }

    @Test
    public void testINSStart_pref() throws Exception{
        int fileFeatStart = 10499052;
        tstStart("chr20", fileFeatStart, VariantContext.Type.INDEL, fileFeatStart - 1 + 1);
    }

    @Test
    public void testINSStart_nopref() throws Exception{
        int fileFeatStart = 10499352;
        tstStart("chr20", fileFeatStart, VariantContext.Type.INDEL, fileFeatStart - 1);
    }

    private VCFVariant tstStart(String chr, int fileFeatStart, VariantContext.Type variantType, int expFeatStart) throws Exception{
        VCFVariant variant = get_hc_mod_Feat(chr, fileFeatStart - 5, fileFeatStart + 5);
        assertEquals(variantType.toString(), variant.getType());
        assertEquals(expFeatStart, variant.getStart());
        return variant;
    }
}
