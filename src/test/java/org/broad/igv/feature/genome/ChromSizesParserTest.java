package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.load.ChromSizesParser;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 4/16/13
 *         Time: 2:43 PM
 */
public class ChromSizesParserTest {

    String testFile = TestUtils.DATA_DIR + "genomes/hg19.chrom.sizes";

    @Test
    public void testParse() throws Exception {

        int expectedCount = 37;

        List<Chromosome> chromosomes = ChromSizesParser.parse(testFile);
        assertEquals(expectedCount, chromosomes.size());

        // chr10	135534747	/gbdb/hg19/hg19.2bit
        Chromosome chr10 = chromosomes.get(9);
        assertEquals("chr10", chr10.getName());
        assertEquals(135534747, chr10.getLength());
    }
}
