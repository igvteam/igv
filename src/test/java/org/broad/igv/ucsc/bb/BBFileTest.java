package org.broad.igv.ucsc.bb;

import htsjdk.tribble.Feature;
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

        BBFile reader = BBFile.openFile(bbFile);

        BBHeader header = reader.readHeader();

        assertNotNull(header);

    }

    @Test
    public void testBigGenePred() throws IOException {

        String path = "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.catLiftOffGenesV1/catLiftOffGenesV1.bb";

        BBFile bbReader = BBFile.openFile(path);
        BBDataSource bbSource = new BBDataSource(bbReader, null);

        String chr = "chr21";
        int start = 26490012;
        int end = 42182827;

        Iterator<Feature> iter = bbSource.getFeatures(chr, start, end);

        int count = 0;
        while (iter.hasNext()) {
            Feature f = iter.next();
            assertEquals(chr, f.getChr());
            assertTrue(f.getStart() <= end && f.getEnd() >= start);
            count++;
        }
        assertTrue(count > 0);

    }


}