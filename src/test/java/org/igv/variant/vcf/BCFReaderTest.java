package org.igv.variant.vcf;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.igv.AbstractHeadlessTest;
import org.igv.exceptions.DataLoadException;
import org.igv.track.FeatureSource;
import org.igv.track.TrackLoader;
import org.igv.track.TribbleFeatureSource;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

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

    /**
     * BCF other than uncompressed version 2.1 fails with an explanation, not a header parsing error (issue #1495).
     * The test files are made from ex2.bcf (version 2.1, uncompressed).
     */
    @Test
    public void unsupportedBCFExplained() throws Exception {
        byte[] v21 = Files.readAllBytes(Paths.get(TestUtils.DATA_DIR + "bcf/ex2.bcf"));
        byte[] v22 = v21.clone();
        v22[4] = 2;     // minor version
        File dir = new File(TestUtils.TMP_OUTPUT_DIR, "bcf");
        dir.mkdirs();

        assertExplained(writeBCF(dir, "v22.bcf", v22, false), "This is a BCF version 2.2 file");
        assertExplained(writeBCF(dir, "v21_compressed.bcf", v21, true),
                "This is a compressed BCF file, but its name does not end in .gz or .bgz");
        assertExplained(writeBCF(dir, "v22_compressed.bcf", v22, true), "This is a compressed BCF version 2.2 file");
        assertEquals(1, new TrackLoader().load(new ResourceLocator(TestUtils.DATA_DIR + "bcf/ex2.bcf"), genome).size());

        // htsjdk decompresses based on the file name, so compressed BCF 2.1 named .gz or .bgz is readable
        assertEquals(1, new TrackLoader().load(new ResourceLocator(writeBCF(dir, "v21_compressed.bcf.gz", v21, true)), genome).size());
        assertEquals(1, new TrackLoader().load(new ResourceLocator(writeBCF(dir, "v21_compressed.bcf.bgz", v21, true)), genome).size());
        assertExplained(writeBCF(dir, "v22_compressed.bcf.gz", v22, true), "This is a compressed BCF version 2.2 file");

        // Without a .bcf extension, the format is recognized from the file contents
        assertExplained(writeBCF(dir, "v22_compressed_no_extension", v22, true), "This is a compressed BCF version 2.2 file");
        assertEquals(1, new TrackLoader().load(new ResourceLocator(writeBCF(dir, "v21_no_extension", v21, false)), genome).size());
    }

    /**
     * The suggested bcftools command is built from the URL path, which must not include a query string such as a
     * signed URL's signature.
     */
    @Test
    public void urlPathExcludesQuery() {
        assertEquals("https://example.org/data/x.bcf",
                new ResourceLocator("https://example.org/data/x.bcf?X-Amz-Signature=abc").getURLPath());
    }

    private static String writeBCF(File dir, String name, byte[] bytes, boolean compressed) throws IOException {
        File file = new File(dir, name);
        try (OutputStream out = compressed ? new BlockCompressedOutputStream(file) : new FileOutputStream(file)) {
            out.write(bytes);
        }
        return file.getAbsolutePath();
    }

    private void assertExplained(String path, String expected) {
        try {
            new TrackLoader().load(new ResourceLocator(path), genome);
            fail("Expected an error loading " + path);
        } catch (DataLoadException e) {
            assertTrue(e.getMessage(), e.getMessage().contains(expected));
            String name = new File(path).getName();
            String vcfName = name.replaceFirst("\\.bcf(\\.gz|\\.bgz)?$", "") + ".vcf.gz";
            assertTrue(e.getMessage(), e.getMessage().contains("bcftools view -Oz -o " + vcfName + " " + name));
        }
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
