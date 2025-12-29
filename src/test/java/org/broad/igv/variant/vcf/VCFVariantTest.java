package org.broad.igv.variant.vcf;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Ignore;
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

    /**
     * Just verifies IGV will read a file marked vcf4.4 in the header, actual format is 4.1.
     *
     * @return
     * @throws Exception
     */
   @Ignore
   @Test
    public void tstFakeV4() throws Exception{

        System.setProperty("samjdk.optimistic_vcf_4_4", "true");  // Doesn't seem to work from command line (gradlew test)

        String filePath = TestUtils.DATA_DIR + "vcf/fake_v4.vcf";
        TestUtils.createIndex(filePath);
        TribbleFeatureSource src = TribbleFeatureSource.getFeatureSource(new ResourceLocator(filePath), genome);

        VCFVariant variant = (VCFVariant) (src.getFeatures("chr20", 10499352 - 5, 10499352 + 5)).next();
        assertEquals(VariantContext.Type.INDEL.toString(), variant.getType());
        assertEquals(10499351, variant.getStart());
    }
}
