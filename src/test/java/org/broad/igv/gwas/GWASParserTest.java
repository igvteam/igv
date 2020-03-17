/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.gwas;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
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
        Map<String, List<GWASFeature>> data = parser.parse();
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
        Map<String, List<GWASFeature>> data = parser.parse();
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
        try {
            Map<String, List<GWASFeature>> data = parser.parse();
        } catch (ParserException e) {
            excepted = true;
        }
        assertTrue(excepted);

    }
}
