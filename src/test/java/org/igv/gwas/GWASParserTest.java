package org.igv.gwas;

import org.igv.AbstractHeadlessTest;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;
import java.util.Map;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012/04/18
 */
public class GWASParserTest extends AbstractHeadlessTest {

    @Test
    public void testParse_underflows() throws Exception {
        GWASParser parser = new GWASParser(new ResourceLocator(TestUtils.DATA_DIR + "gwas/smallp.gwas"), genome);
        GWASData data = parser.parse();
        for (List<GWASFeature> features : data.values()) {
            for (GWASFeature f : features) {
                double val = f.value;
                assertFalse("Value is infinite", Double.isInfinite(val));
                assertFalse("Value isNan", Double.isNaN(val));
                assertFalse("Value is 0", val == 0.0f);
            }
        }
    }

    @Test
    public void testParseStarts() throws Exception {
        GWASParser parser = new GWASParser(new ResourceLocator(TestUtils.DATA_DIR + "gwas/smallp.gwas"), genome);
        GWASData data = parser.parse();
        List<GWASFeature> features = data.get("chr6");

        int[] expStarts = {29622220, 29623739, 29623739};
        for (int ii = 0; ii < expStarts.length; ii++) {
            assertEquals(expStarts[ii], features.get(ii).position);
        }
    }

    @Test
    public void testBadPs() throws Exception {
        String[] finames = {"badp_neg.gwas", "badp_text.gwas", "badp_zero.gwas"};
        for (String finame : finames) {
            tstParseBad(finame);
        }
    }

    public void tstParseBad(String finame) throws Exception {
        GWASParser parser = new GWASParser(new ResourceLocator(TestUtils.DATA_DIR + "gwas/" + finame), genome);
        boolean excepted = false;

        GWASData data = parser.parse();

        assertTrue(data.isEmpty());

    }
}
