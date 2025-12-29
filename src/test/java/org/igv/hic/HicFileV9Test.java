package org.igv.hic;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class HicFileV9Test {

    static String testURL = "https://www.encodeproject.org/files/ENCFF053VBX/@@download/ENCFF053VBX.hic";
    static HicFile hicFile;

    @BeforeClass
    public static void setUp() throws Exception {
        hicFile = new HicFile(testURL, null);
    }

    @Test
    public void getVersion() {
        int version = hicFile.getVersion();
        assertEquals(9, version);
    }


    @Test
    public void getNormalizationTypes() throws IOException {
        List<String> normTypes = hicFile.getNormalizationTypes();
        assertEquals(6, normTypes.size());
    }


    /**
     * Known HIC v9 files:
     * https://www.encodeproject.org/files/ENCFF053VBX/@@download/ENCFF053VBX.hic is version 9
     * https://www.encodeproject.org/files/ENCFF876OWE/@@download/ENCFF876OWE.hic is version 9
     * https://www.encodeproject.org/files/ENCFF685BLG/@@download/ENCFF685BLG.hic is version 9
     * https://www.encodeproject.org/files/ENCFF778OYA/@@download/ENCFF778OYA.hic is version 9
     * https://www.encodeproject.org/files/ENCFF020DPP/@@download/ENCFF020DPP.hic is version 9
     * https://www.encodeproject.org/files/ENCFF216QQM/@@download/ENCFF216QQM.hic is version 9
     * https://www.encodeproject.org/files/ENCFF542BHD/@@download/ENCFF542BHD.hic is version 9
     * https://www.encodeproject.org/files/ENCFF750AOC/@@download/ENCFF750AOC.hic is version 9
     * https://www.encodeproject.org/files/ENCFF301BWY/@@download/ENCFF301BWY.hic is version 9
     * https://www.encodeproject.org/files/ENCFF965PEE/@@download/ENCFF965PEE.hic is version 9
     * https://www.encodeproject.org/files/ENCFF174LAF/@@download/ENCFF174LAF.hic is version 9
     * https://www.encodeproject.org/files/ENCFF234MZQ/@@download/ENCFF234MZQ.hic is version 9
     * https://www.encodeproject.org/files/ENCFF080DPJ/@@download/ENCFF080DPJ.hic is version 9
     */
}
