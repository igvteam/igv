package org.broad.igv.data;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Aug-12
 */
public class WiggleParserTest extends AbstractHeadlessTest{

    @Test
    public void testParseBedgraph_onefeat() throws Exception {
        String filepath = TestUtils.DATA_DIR + "wig/test_coords.bedgraph";
        WiggleParser parser = new WiggleParser(new ResourceLocator(filepath));
        String chr = "chr1";
        WiggleDataset ds = parser.parse();
        assertArrayEquals(new String[]{chr}, ds.getChromosomes());
        assertEquals(1.0, ds.getDataMax(), 0.01);
        assertArrayEquals(new int[]{1}, ds.getStartLocations(chr));
    }

    @Test
    public void testParseBedgraph_stats() throws Exception {

        String filepath = TestUtils.DATA_DIR + "wig/test_stats.bedgraph";
        WiggleParser parser = new WiggleParser(new ResourceLocator(filepath));

        String chr = "chr1";
        WiggleDataset ds = parser.parse();

        int len = 998;
        assertEquals(len, ds.getStartLocations(chr).length);
        int[] starts = new int[len];
        int[] ends = new int[len];
        Arrays.fill(starts, 0);
        Arrays.fill(ends, 1000);
        float[] vals = new float[]{0.15244511f, 0.9837612f, 0.17782128f, 0.81641054f, 0.89790136f};
        int[] fiStarts = ds.getStartLocations(chr);
        int[] fiEnds = ds.getEndLocations(chr);
        float[] fiVals = ds.getData("", chr);
        for(int ii=0; ii < len; ii++){
            assertEquals(starts[ii], fiStarts[ii]);
            assertEquals(ends[ii], fiEnds[ii]);
            if(ii < vals.length){
                assertEquals(vals[ii], fiVals[ii], 0.0000001f);
            }
        }

        assertEquals(0.9998707, ds.dataMax, 0.000001);

        //expected values calculated with excel which uses a different algorithm
        //Since we actually interpolate/estimate, don't expect too much precision
        assertEquals(0.086, ds.percent10, 0.001);
        assertEquals(0.89, ds.percent90, 0.01);
    }


}
