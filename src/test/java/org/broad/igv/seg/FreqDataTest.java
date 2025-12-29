package org.broad.igv.seg;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;

/**
 * @author jrobinso
 * @date Oct 13, 2010
 */
public class FreqDataTest extends AbstractHeadlessTest {

    @Test
    public void test() throws IOException {

        Genome genome = TestUtils.loadGenome();

        String segfile = TestUtils.DATA_DIR + "seg/canFam2_hg18.seg";

        SegmentedDataSet sd = SegmentFileParser.loadSegments(new ResourceLocator(segfile), genome);

        FreqData fd = new FreqData(sd, genome);

        Assert.assertNotNull(fd);
    }
}
