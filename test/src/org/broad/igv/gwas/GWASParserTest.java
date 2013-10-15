/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
