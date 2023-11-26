package org.broad.igv.feature.genome;

import org.broad.igv.feature.Cytoband;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

public class CytobandSourceBBTest {

    @Test
    public void getCytobands() throws IOException {
        String url = TestUtils.DATA_DIR + "bb/cytoBandMapped.bb";
        Genome genome = null;
        CytobandSourceBB src = new CytobandSourceBB(url, null);
        List<Cytoband> cytobands = src.getCytobands("chr1");
        Cytoband last = cytobands.get(cytobands.size() - 1);
        assertEquals(248387328, last.getEnd());
    }
}