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

package org.broad.igv.variant.vcf;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
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
