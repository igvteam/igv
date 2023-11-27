package org.broad.igv.ucsc;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.Map;

import static org.junit.Assert.*;

public class TrixTest {
    @Test
    public void testTrix() throws IOException {
        String ixFile = TestUtils.DATA_DIR + "bb/ixIxx/GCF_000009045.1_ASM904v1.ncbiGene.ix";
        String ixxFile = TestUtils.DATA_DIR + "bb/ixIxx/GCF_000009045.1_ASM904v1.ncbiGene.ixx";
        Trix trix = new Trix(ixxFile, ixFile);
        Map<String, String[]> results = trix.search("ykoX");
        String[] exactMatches = results.get("ykox");
        assertEquals("NP_389226.1", exactMatches[0]);
    }
}