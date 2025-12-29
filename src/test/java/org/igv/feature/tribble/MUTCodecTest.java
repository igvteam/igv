package org.igv.feature.tribble;

import org.igv.AbstractHeadlessTest;
import org.igv.feature.Mutation;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 * Date: 4/8/13
 * Time: 2:51 PM
 */
public class MUTCodecTest extends AbstractHeadlessTest {


    @Test
    public void testSampleMeta() throws Exception {

        String path = TestUtils.DATA_DIR + "maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf";
        assertTrue(MUTCodec.isMutationAnnotationFile(new ResourceLocator(path)));

        MUTCodec codec = new MUTCodec(path, genome);
        int expectedSampleCount = 146;
        String[] samples = codec.getSamples();
        assertEquals(expectedSampleCount, samples.length);
        assertTrue(samples[0].startsWith("TCGA"));
    }

    @Test
    public void testBroadMaf() throws Exception {
        String path = TestUtils.DATA_DIR + "maf/test.maf";
        assertTrue(MUTCodec.isMutationAnnotationFile(new ResourceLocator(path)));

        MUTCodec codec = new MUTCodec(path, genome);

        BufferedReader bufferedReader = null;
        try {
            bufferedReader = new BufferedReader(new FileReader(path));
            String nextLine;
            ArrayList<Mutation> mutations = new ArrayList<>();
            while ((nextLine = bufferedReader.readLine()) != null) {

                Mutation mut = codec.decode(nextLine);
                if (mut != null) {
                    mutations.add(mut);
                }
            }
            assertEquals(29, mutations.size());
        } finally {
            bufferedReader.close();
        }
    }

    @Test
    public void testTcgaMaf() throws Exception {
        String path = TestUtils.DATA_DIR + "maf/tcga_test.maf";
        assertTrue(MUTCodec.isMutationAnnotationFile(new ResourceLocator(path)));

        MUTCodec codec = new MUTCodec(path, genome);

        BufferedReader bufferedReader = null;
        try {
            bufferedReader = new BufferedReader(new FileReader(path));
            String nextLine;
            ArrayList<Mutation> mutations = new ArrayList<>();
            while ((nextLine = bufferedReader.readLine()) != null) {

                Mutation mut = codec.decode(nextLine);
                if (mut != null) {
                    mutations.add(mut);
                }
            }
            assertEquals(17, mutations.size());
        } finally {
            bufferedReader.close();
        }
    }

    @Test
    public void testMafLite() throws Exception {
        String path = TestUtils.DATA_DIR + "maf/test.maflite.tsv";
        assertTrue(MUTCodec.isMutationAnnotationFile(new ResourceLocator(path)));

        MUTCodec codec = new MUTCodec(path, genome);

        BufferedReader bufferedReader = null;
        try {
            bufferedReader = new BufferedReader(new FileReader(path));
            String nextLine;
            ArrayList<Mutation> mutations = new ArrayList<>();
            while ((nextLine = bufferedReader.readLine()) != null) {

                Mutation mut = codec.decode(nextLine);
                if (mut != null) {
                    mutations.add(mut);
                }
            }
            assertEquals(29, mutations.size());
        } finally {
            bufferedReader.close();
        }
    }

}
