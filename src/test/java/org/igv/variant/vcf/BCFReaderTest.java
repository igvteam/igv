package org.igv.variant.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedInputStream;
import org.igv.AbstractHeadlessTest;
import org.igv.track.FeatureSource;
import org.igv.track.TrackLoader;
import org.igv.track.TribbleFeatureSource;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * @author jacob
 * @date 2013-Jun-14
 */
public class BCFReaderTest extends AbstractHeadlessTest {


    /**
     * Just test that we load a BCF file without crashing
     *
     * @throws Exception
     */
    @Test
    public void loadBCF() throws Exception {
        String path = TestUtils.DATA_DIR + "bcf/ex2.bcf";
        FeatureSource source = TribbleFeatureSource.getFeatureSource(new ResourceLocator(path), genome);
        Iterator<Feature> features = source.getFeatures("20", 14000, 1300000);
        int count = 0;

        while (features.hasNext()) {
            features.next();
            count++;
        }
        assertTrue("No features read", count > 0);
    }

    @Test
    public void bgzfBcfIsQueriedThroughCsiIndex() throws Exception {
        File bcf = writeIndexedBCF("indexed.bcf");

        TribbleFeatureSource source = TribbleFeatureSource.getFeatureSource(new ResourceLocator(bcf.getPath()), genome);

        assertTrue(source.isIndexed());
        assertEquals(List.of("v2"), ids(source.getFeatures("chr1", 1500, 2500)));
    }

    @Test
    public void bgzfBcfNamedGzIsQueriedThroughCsiIndex() throws Exception {
        File bcf = writeIndexedBCF("renamed.bcf");
        File gz = new File(bcf.getPath() + ".gz");
        assertTrue(bcf.renameTo(gz));
        assertTrue(new File(bcf.getPath() + ".csi").renameTo(new File(gz.getPath() + ".csi")));

        TribbleFeatureSource source = TribbleFeatureSource.getFeatureSource(new ResourceLocator(gz.getPath()), genome);

        assertTrue(source.isIndexed());
        assertEquals(List.of("v2"), ids(source.getFeatures("chr1", 1500, 2500)));
    }

    @Test
    public void bgzfBcfWithoutIndexIsReadThrough() throws Exception {
        File bcf = writeIndexedBCF("unindexed.bcf");
        assertTrue(new File(bcf.getPath() + ".csi").delete());

        TribbleFeatureSource source = TribbleFeatureSource.getFeatureSource(new ResourceLocator(bcf.getPath()), genome);

        assertFalse(source.isIndexed());
        assertEquals(List.of("v2"), ids(source.getFeatures("chr1", 1500, 2500)));
    }

    @Test
    public void bgzfBcfLoadsAsTrack() throws Exception {
        File bcf = writeIndexedBCF("track.bcf");
        assertEquals(1, new TrackLoader().load(new ResourceLocator(bcf.getPath()), genome).size());
    }

    @Test
    public void bgzfBcfWithoutExtensionIsRecognizedByContent() throws Exception {
        File bcf = writeIndexedBCF("no_extension.bcf");
        File noExtension = new File(bcf.getParentFile(), "no_extension");
        assertTrue(bcf.renameTo(noExtension));

        assertEquals(1, new TrackLoader().load(new ResourceLocator(noExtension.getPath()), genome).size());
    }

    /**
     * Writes BCF as bcftools writes it, BGZF-compressed BCF 2.2 with a CSI index, holding SNPs v1, v2 and v3 on chr1
     * at 1000, 2000 and 3000.
     */
    private static File writeIndexedBCF(String name) throws IOException {
        File dir = new File(TestUtils.TMP_OUTPUT_DIR, "bcf");
        dir.mkdirs();
        File bcf = new File(dir, name);

        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 100000)));
        VCFHeader header = new VCFHeader();
        header.setSequenceDictionary(dictionary);
        try (VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputPath(bcf.toPath())
                .setReferenceDictionary(dictionary)
                .setOption(Options.INDEX_ON_THE_FLY)
                .build()) {
            writer.writeHeader(header);
            for (int i = 1; i <= 3; i++) {
                writer.add(new VariantContextBuilder("test", "chr1", i * 1000, i * 1000,
                        List.of(Allele.create("A", true), Allele.create("C"))).id("v" + i).make());
            }
        }

        try (BlockCompressedInputStream in = new BlockCompressedInputStream(bcf.toPath())) {
            assertArrayEquals("BCF\2\2".getBytes(), in.readNBytes(5));
        }
        assertTrue(new File(bcf.getPath() + ".csi").exists());
        return bcf;
    }

    private static List<String> ids(Iterator<Feature> features) {
        List<String> ids = new ArrayList<>();
        features.forEachRemaining(f -> ids.add(((VCFVariant) f).getID()));
        return ids;
    }

    @Test
    public void loadBCFChromAlias() throws Exception {
        String path = TestUtils.DATA_DIR + "bcf/ex2.bcf";
        FeatureSource source = TribbleFeatureSource.getFeatureSource(new ResourceLocator(path), genome);
        Iterator<Feature> features = source.getFeatures("chr20", 14000, 1300000);
        int count = 0;

        while (features.hasNext()) {
            features.next();
            count++;
        }
        assertTrue("No features read", count > 0);
    }



    /**
     * Compare a BCF and VCF file
     *
     * @throws Exception
     */
    @Test
    public void compareBCFtoVCF() throws Exception {
        String BCF2path = TestUtils.DATA_DIR + "bcf/ex2.bcf";
        FeatureSource BCF2source = TribbleFeatureSource.getFeatureSource(new ResourceLocator(BCF2path), genome);
        Iterator<Feature> BCF2features = BCF2source.getFeatures("chr20", 14000, 1300000);
        List<VCFVariant> BCF2List = new ArrayList<VCFVariant>();

        String VCFpath = TestUtils.DATA_DIR + "vcf/ex2.vcf";
        TestUtils.createIndex(VCFpath);
        FeatureSource VCFsource = TribbleFeatureSource.getFeatureSource(new ResourceLocator(VCFpath), genome);
        Iterator<Feature> VCFfeatures = VCFsource.getFeatures("chr20", 14000, 1300000);
        List<VCFVariant> VCFList = new ArrayList<VCFVariant>();

        while (BCF2features.hasNext()) {
            VCFVariant bcfV = (VCFVariant) BCF2features.next();
            VCFVariant vcfV = (VCFVariant) VCFfeatures.next();

            assertEquals(vcfV.getType(), bcfV.getType());

            BCF2List.add(bcfV);
            VCFList.add(vcfV);
        }

        TestUtils.assertFeatureListsEqual(VCFList.iterator(), BCF2List.iterator());
    }

    //Quick method for checking if a bcf file has the magic header
    @Ignore
    //@Test
    public void rawTestFile() throws Exception {
        String path = "/path/to/myfile.bcf";
        PositionalBufferedStream ps = new PositionalBufferedStream(new FileInputStream(path));

        BCF2Codec codec = new BCF2Codec();
        codec.readHeader(ps);

    }


}
