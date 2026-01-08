package org.igv.hic;

import org.igv.feature.genome.Genome;
import org.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class HicFileLocalTest {

    static String dataFile = TestUtils.DATA_DIR + "hic/test_chr22.hic";
    static HicFile hicFile;

    @BeforeClass
    public static void setUp() throws Exception {
        Genome genome = TestUtils.mockUCSCGenome();
        hicFile = new HicFile(dataFile, genome);
    }


    @Test
    public void getVersion() {
        int version = hicFile.getVersion();
        assertEquals(8, version);
    }

    @Test
    public void getNormalizationVector() throws IOException {

        String type = "KR";
        String chr = "22";
        String unit = "BP";
        int binSize = 100000;
        NormalizationVector normVector = hicFile.getNormalizationVector(type, chr, unit, binSize);
        assertEquals(515, normVector.getnValues());

        double[] values =  normVector.getValues(200, 300);
        assertEquals(0.9477306190631148, values[0], 0.00001);

        values = normVector.getValues(0, 100);
        assertTrue(Double.isNaN(values[0]));
    }

}
