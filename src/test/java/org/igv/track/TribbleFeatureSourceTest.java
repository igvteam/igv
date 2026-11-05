package org.igv.track;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.tabix.TabixIndexType;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.igv.AbstractHeadlessTest;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.igv.variant.Variant;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class TribbleFeatureSourceTest extends AbstractHeadlessTest {

    @Test
    public void bgzippedVcfWithOnlyCsiIndexIsQueriedThroughIt() throws Exception {
        File vcf = writeIndexedVcf("csi_only.vcf.gz", TabixIndexType.CSI);
        assertTrue(new File(vcf.getPath() + ".csi").exists());
        assertFalse(new File(vcf.getPath() + ".tbi").exists());

        TribbleFeatureSource source = TribbleFeatureSource.getFeatureSource(new ResourceLocator(vcf.getPath()), genome);

        assertTrue(source.isIndexed());
        assertEquals(List.of("v2"), ids(source.getFeatures("chr1", 1500, 2500)));
    }

    @Test
    public void bgzippedVcfWithTbiIndexIsQueriedThroughIt() throws Exception {
        File vcf = writeIndexedVcf("tbi_only.vcf.gz", TabixIndexType.TBI);

        TribbleFeatureSource source = TribbleFeatureSource.getFeatureSource(new ResourceLocator(vcf.getPath()), genome);

        assertTrue(source.isIndexed());
        assertEquals(List.of("v2"), ids(source.getFeatures("chr1", 1500, 2500)));
    }

    /**
     * Writes a bgzipped VCF with SNPs v1, v2 and v3 on chr1 at 1000, 2000 and 3000, indexed while writing.
     */
    private static File writeIndexedVcf(String name, TabixIndexType indexType) {
        File dir = new File(TestUtils.TMP_OUTPUT_DIR, "tribble_feature_source");
        dir.mkdirs();
        File vcf = new File(dir, name);

        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 100000)));
        VCFHeader header = new VCFHeader();
        header.setSequenceDictionary(dictionary);

        try (VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputPath(vcf.toPath())
                .setReferenceDictionary(dictionary)
                .setOption(Options.INDEX_ON_THE_FLY)
                .setTabixIndexType(indexType)
                .build()) {
            writer.writeHeader(header);
            for (int i = 1; i <= 3; i++) {
                writer.add(new VariantContextBuilder("test", "chr1", i * 1000, i * 1000,
                        List.of(Allele.create("A", true), Allele.create("C"))).id("v" + i).make());
            }
        }
        return vcf;
    }

    private static List<String> ids(Iterator<Feature> features) {
        List<String> ids = new ArrayList<>();
        features.forEachRemaining(f -> ids.add(((Variant) f).getID()));
        return ids;
    }
}
