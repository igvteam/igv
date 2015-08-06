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
import org.broad.igv.util.collections.DoubleArrayList;
import org.broad.igv.util.collections.IntArrayList;
import org.junit.Test;

import java.util.LinkedHashMap;
import java.util.List;

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
        LinkedHashMap<String, DoubleArrayList> values = data.getValues();
        for (String chr : values.keySet()) {
            DoubleArrayList floats = values.get(chr);
            for (int ff = 0; ff < floats.size(); ff++) {
                double val = floats.get(ff);
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
        IntArrayList startLocs = data.getLocations().get("chr6");

        int[] expStarts = {29622220,29623739,29623739};
        for(int ii=0; ii < expStarts.length; ii++){
            assertEquals(expStarts[ii], startLocs.get(ii));
        }

    }

    @Test
    public void testBadPs() throws Exception {
        String[] finames = {"badp_neg.gwas", "badp_text.gwas", "badp_zero.gwas"};
        for (String finame : finames) {
            tstParseBad(finame);
        }
    }

    @Test
    public void testUnsorted() throws Exception{
        String path = "random.gwas";
        tstParseBad(path);
    }

    public void tstParseBad(String finame) throws Exception {
        GWASParser parser = new GWASParser(new ResourceLocator(TestUtils.DATA_DIR + "gwas/" + finame), genome);
        boolean excepted = false;
        try {
            GWASData data = parser.parse();
        } catch (ParserException e) {
            excepted = true;
        }
        assertTrue(excepted);

    }

    @Test
    public void testLoadGWAS() throws Exception{
        ResourceLocator locator = new ResourceLocator(TestUtils.DATA_DIR + "gwas/smallp.gwas");
        List<Track> tracks = (new TrackLoader()).load(locator, genome);
        GWASTrack track = (GWASTrack) tracks.get(0);
        String desc = track.getDescription("chr6", 1);

        String[] lines = desc.split("<br>");
        String[] expTokens = new String[]{"rs29228", "6", "29623739", "0.931148124684"};
        int offset = 3;
        for(int tn=0; tn < expTokens.length; tn++){
            String token = lines[tn+offset];
            String[] areas = token.split("\\s");
            String value = areas[1];
            assertEquals("Value for field " + areas[0] + " not equal", expTokens[tn], value);
        }
    }
}
