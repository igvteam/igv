package org.broad.tribble;


import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * User: jrobinso
 * Date: Jul 17, 2010
 * Time: 9:28:49 PM
 */
@Ignore("Requires largedata bundle")
public class TribbleFeatureReaderTest {

    static String testFile = TestUtils.LARGE_DATA_DIR + "CEU.SRP000032.2010_03_v4.0.genotypes.head.vcf";
    static AbstractFeatureReader<VariantContext, ?> bfr;

    @BeforeClass
    public static void setUpClass() throws IOException {
        File idxFile = new File(testFile + ".idx");
        FeatureCodec codec = new VCFCodec();
        createIndex(idxFile, codec);

        bfr = AbstractFeatureReader.getFeatureReader(testFile, codec);
    }

    @AfterClass
    public static void tearDownClass() throws IOException {
        bfr.close();
    }

    @Test @Ignore("Requires largedata bundle")
    public void testQuery() throws IOException {

        String chr = "1";
        int start = 1718546;
        int end = 1748915;
        int[] expectedStarts = {1718547, 1718829, 1723079, 1724830, 1731376, 1733967, 1735586, 1736016, 1738594,
                1739272, 1741124, 1742815, 1743224, 1748886, 1748914};


        Iterator<VariantContext> iter = bfr.query(chr, start, end);
        int count = 0;
        while (iter.hasNext()) {
            VariantContext feat = iter.next();
            int expStart = expectedStarts[count];
            assertEquals(expStart, feat.getStart());
            count++;
        }

        Assert.assertEquals(15, count);

    }

    @Test @Ignore("Requires largedata bundle")
    public void testGetSequenceNames() throws Exception {
        Set<String> expectedSequences = new HashSet(Arrays.asList("1", "2"));

        int count = 0;
        for (String s : bfr.getSequenceNames()) {
            assertTrue(expectedSequences.contains(s));
            count++;
        }

        assertEquals(expectedSequences.size(), count);
    }


    private static void createIndex(File idxFile, FeatureCodec codec) throws IOException {
        if (idxFile.exists()) {
            idxFile.delete();
        }

        // Create the index  
        Index idx = IndexFactory.createIntervalIndex(new File(testFile), codec, 10);

        LittleEndianOutputStream stream = null;
        try {
            stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(idxFile)));
            idx.write(stream);
        } finally {
            if (stream != null) {
                stream.close();
            }
        }
        idxFile.deleteOnExit();

    }

}