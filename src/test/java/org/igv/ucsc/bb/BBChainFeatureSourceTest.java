package org.igv.ucsc.bb;

import org.igv.feature.Exon;
import org.igv.feature.IGVFeature;
import org.igv.feature.PSLRecord;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.*;

import static org.junit.Assert.*;

public class BBChainFeatureSourceTest {

    // Example files from https://genome.ucsc.edu/goldenPath/help/bigChain.html -- 4 chains, 11 links
    @Test
    public void testBigChain() throws IOException {

        BBFile chainReader = new BBFile(TestUtils.DATA_DIR + "bb/bigChain.bb", null);
        BBFile linkReader = new BBFile(TestUtils.DATA_DIR + "bb/bigChain.link.bb", null);
        BBChainFeatureSource source = new BBChainFeatureSource(chainReader, linkReader, null);

        String chr = "chr22_KI270731v1_random";
        Iterator<IGVFeature> iter = source.getFeatures(chr, 0, Integer.MAX_VALUE);

        Map<String, PSLRecord> chains = new HashMap<>();
        int blockCount = 0;
        while (iter.hasNext()) {
            IGVFeature f = iter.next();
            assertTrue(f instanceof PSLRecord);
            PSLRecord chain = (PSLRecord) f;
            chains.put(chain.getIdentifier(), chain);
            for (Exon exon : chain.getExons()) {
                assertTrue(exon.getStart() >= chain.getStart() && exon.getEnd() <= chain.getEnd());
                blockCount++;
            }
        }
        assertEquals(4, chains.size());
        assertEquals(11, blockCount);

        // chain 569  chr22_KI270731v1_random:28198-62665  qName=chr3 qSize=160039680 qStart=101389330 qEnd=101422472
        PSLRecord chain = chains.get("569");
        assertEquals("chr3", chain.getName());
        assertEquals(28198, chain.getStart());
        assertEquals(62665, chain.getEnd());
        assertEquals(150754, chain.getTSize());
        assertEquals(160039680, chain.getqSize());
        assertTrue(chain.getExons().size() > 1);
        assertEquals(28198, chain.getExons().get(0).getStart());
        assertEquals(28199, chain.getExons().get(0).getEnd());

        // Only links overlapping the query interval are attached
        iter = source.getFeatures(chr, 40000, 41000);
        chain = (PSLRecord) iter.next();
        assertEquals("569", chain.getIdentifier());
        for (Exon exon : chain.getExons()) {
            assertTrue(exon.getEnd() > 40000 && exon.getStart() < 41000);
        }
        assertFalse(iter.hasNext());

        // A chain overlapping an interval that contains none of its links gets 1 bp placeholders at the interval edges
        iter = source.getFeatures(chr, 30000, 31000);
        chain = (PSLRecord) iter.next();
        assertEquals("569", chain.getIdentifier());
        assertEquals(2, chain.getExons().size());
        assertEquals(30000, chain.getExons().get(0).getStart());
        assertEquals(30999, chain.getExons().get(1).getStart());
    }

    @Test
    public void testNoLinkFile() throws IOException {
        BBFile chainReader = new BBFile(TestUtils.DATA_DIR + "bb/bigChain.bb", null);
        BBChainFeatureSource source = new BBChainFeatureSource(chainReader, null, null);
        Iterator<IGVFeature> iter = source.getFeatures("chr22_KI270731v1_random", 0, Integer.MAX_VALUE);
        int count = 0;
        while (iter.hasNext()) {
            PSLRecord chain = (PSLRecord) iter.next();
            assertFalse(chain.hasExons());
            count++;
        }
        assertEquals(4, count);
    }
}
