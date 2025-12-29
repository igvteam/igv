package org.broad.igv.variant.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

import static junit.framework.Assert.assertEquals;


/**
 * @author jacob
 * @date 2012/05/11
 */
public class VCFWriterTest extends AbstractHeadlessTest {

    String inpath = TestUtils.DATA_DIR + "vcf/SRP32_v4.sorted.0.vcf";
    File outFile = new File(TestUtils.TMP_OUTPUT_DIR, "testwriterout.vcf");

    private VariantContextWriter getWriter() {
        SAMSequenceDictionary seqDict = new SAMSequenceDictionary();
        EnumSet<Options> options = VariantContextWriterBuilder.DEFAULT_OPTIONS;
        options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(seqDict)
                .setOptions(options);

        VariantContextWriter writer = builder.setOutputFile(outFile).build();
        return writer;
    }

    @Test
    public void testWriteHeader() throws Exception {

        FeatureCodec codec = CodecFactory.getCodec(inpath, genome);
        AbstractFeatureReader<Feature, ?> bfs = AbstractFeatureReader.getFeatureReader(inpath, codec, false);
        VCFHeader header0 = (VCFHeader) bfs.getHeader();

        VariantContextWriter writer = getWriter();
        writer.writeHeader(header0);
        writer.close();

        AbstractFeatureReader<Feature, ?> bfs1 = AbstractFeatureReader.getFeatureReader(outFile.getAbsolutePath(), codec, false);
        VCFHeader header1 = (VCFHeader) bfs1.getHeader();

        assertHeadersEquals(header0, header1);
    }

    @Test
    public void testWriteRecords() throws Exception {

        FeatureCodec codec = CodecFactory.getCodec(inpath, genome);
        AbstractFeatureReader<VCFVariant, ?> bfs = AbstractFeatureReader.getFeatureReader(inpath, codec, false);
        Iterable<VCFVariant> iter0 = bfs.iterator();
        List<VCFVariant> list0 = new ArrayList<VCFVariant>();
        VCFHeader header0 = (VCFHeader) bfs.getHeader();

        VariantContextWriter writer = getWriter();
        writer.writeHeader(header0);
        for (VCFVariant var : iter0) {
            writer.add(var.getVariantContext());
            list0.add(var);
        }
        writer.close();

        AbstractFeatureReader<VCFVariant, ?> bfs1 = AbstractFeatureReader.getFeatureReader(outFile.getAbsolutePath(), codec, false);
        Iterable<VCFVariant> iter1 = bfs1.iterator();
        VCFHeader header1 = (VCFHeader) bfs1.getHeader();

        assertHeadersEquals(header0, header1);
        int n = 0;

        for (VCFVariant var1 : iter1) {
            VCFVariant var0 = list0.get(n++);

            assertVCFVariantsEqual(var0, var1);
        }

    }

    public static void assertVCFVariantsEqual(VCFVariant var0, VCFVariant var1){
        TestUtils.assertFeaturesEqual(var0, var1);

        assertEquals(var0.getType(), var1.getType());
        assertEquals(var0.getID(), var1.getID());
        assertEquals(var0.getSampleNames(), var1.getSampleNames());
        assertEquals(var0.getAttributes(), var1.getAttributes());
    }

    public static void assertHeadersEquals(VCFHeader header0, VCFHeader header1) {


        assertEquals(header0.getColumnCount(), header1.getColumnCount());
        assertEquals(header0.getGenotypeSamples(), header1.getGenotypeSamples());
        assertEquals(header0.getContigLines(), header1.getContigLines());

        for (VCFInfoHeaderLine line0 : header0.getInfoHeaderLines()) {
            VCFInfoHeaderLine line1 = header1.getInfoHeaderLine(line0.getID());
            assertEquals(line0.getCount(), line1.getCount());
            assertEquals(line0.getType(), line1.getType());
            assertEquals(line0.getDescription(), line1.getDescription());
            assertEquals(0, line0.compareTo(line1));
        }
    }
}
