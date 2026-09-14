package org.igv.variant.vcf;

import htsjdk.tribble.Feature;
import org.igv.AbstractHeadlessTest;
import org.igv.track.TrackLoader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.igv.variant.Variant;
import org.igv.variant.VariantTrack;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * gVCF records: the symbolic <NON_REF> and <*> alleles stand for "any other allele", not sequence.
 */
public class GVCFTest extends AbstractHeadlessTest {

    private static final String GVCF = TestUtils.DATA_DIR + "vcf/gvcf_nonref.gvcf";

    /**
     * A gVCF indel covers the same bases as the same indel in a VCF -- the padding base shared by the reference and
     * the real alternate allele is skipped, whatever the <NON_REF> or <*> allele (issue #1620).  Coordinates are
     * 0-based, end exclusive.
     */
    @Test
    public void testIndelCoordinates() {
        List<Variant> variants = variants(GVCF);

        assertCoordinates(variants.get(1), 1007, 1008);     // AT -> A,<NON_REF>: the deleted T
        assertCoordinates(variants.get(3), 1020, 1022);     // ACG -> A,<NON_REF>: the deleted CG
        assertCoordinates(variants.get(4), 1030, 1031);     // AT -> A,<*>
        assertCoordinates(variants.get(5), 1040, 1040);     // A -> ATT,<NON_REF>: an insertion after the A
    }

    /**
     * SNPs and reference blocks are unchanged.
     */
    @Test
    public void testSnpAndReferenceBlock() {
        List<Variant> variants = variants(GVCF);
        assertCoordinates(variants.get(0), 999, 1006);      // A -> <NON_REF>, END=1006
        assertCoordinates(variants.get(2), 1009, 1010);     // C -> G,<NON_REF>
    }

    /**
     * A bgzipped gVCF with a ".gvcf.gz" extension loads (issue #1388).
     */
    @Test
    public void testGvcfGzExtension() {
        String path = GVCF + ".gz";
        assertEquals("gvcf", new ResourceLocator(path).getFormat());
        assertEquals(6, variants(path).size());
    }

    private void assertCoordinates(Variant variant, int start, int end) {
        String label = variant.getReference() + " -> " + variant.getAlternateAlleles();
        assertEquals(label + " start", start, variant.getStart());
        assertEquals(label + " end", end, variant.getEnd());
    }

    private List<Variant> variants(String path) {
        VariantTrack track = (VariantTrack) new TrackLoader().load(new ResourceLocator(path), genome).get(0);
        List<Feature> features = track.getFeatures("chr1", 0, 5000);
        return features.stream().map(f -> (Variant) f).toList();
    }
}
