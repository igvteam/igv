package org.broad.igv.ucsc.bb;

import htsjdk.tribble.Feature;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.Iterator;

import static org.junit.Assert.*;
import static org.junit.Assert.assertTrue;

public class BBFileTest {

    @Test
    public void testGene() throws IOException {

        String bbFile = TestUtils.DATA_DIR + "bb/GCF_000009045.1_ASM904v1.ncbiGene.bb";
        BBFile reader = new BBFile(bbFile, null);
        BBHeader header = reader.readHeader();
        assertNotNull(header);
    }


    @Test
    public void testExtraIndexSearch() throws IOException {

        String path = TestUtils.DATA_DIR + "genomes/GCF_000002655.1.chromAlias.bb";
        BBFile bbReader = new BBFile(path, null);

        // There are 5 extra indexes, 1 for each alias
        String ncbiName = "3";
        BasicFeature f1 = (BasicFeature) bbReader.search(ncbiName);
        assertNotNull(f1);
        assertEquals(ncbiName, f1.getAttribute("ncbi"));

        String ucscName = "chr2";
        BasicFeature f2 = (BasicFeature) bbReader.search(ucscName);
        assertEquals(ucscName, f2.getAttribute("ucsc"));
    }


}