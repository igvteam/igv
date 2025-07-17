package org.broad.igv.feature.genome;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import static org.junit.Assert.*;

public class GenomeDownloadUtilsTest {

    @Test
    public void isAnnotationsDownloadable() {
        assertTrue(GenomeDownloadUtils.isAnnotationsDownloadable(TestUtils.DATA_DIR +  "genomes/json/hg38_twobit.json"));

    }

    @Test
    public void isSequenceDownloadable() {

        assertTrue(GenomeDownloadUtils.isSequenceDownloadable(TestUtils.DATA_DIR +  "genomes/json/hg38_twobit.json"));
        //assertFalse(GenomeDownloadUtils.isSequenceDownloadable(TestUtils.DATA_DIR +  "genomes/json/hg18.unittest.json"));
        //assertFalse(GenomeDownloadUtils.isSequenceDownloadable(TestUtils.DATA_DIR +  "genomes/json/sacCer3.json"));

    }
}