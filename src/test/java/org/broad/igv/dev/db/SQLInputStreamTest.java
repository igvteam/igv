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

package org.broad.igv.dev.db;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.example.ExampleBinaryCodec;
import htsjdk.tribble.readers.PositionalBufferedStream;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.sql.ResultSet;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Aug-22
 */
public class SQLInputStreamTest extends AbstractHeadlessTest {

    @Test
    public void basicTestReadBinary() throws Exception {
        ResultSet rs = DBManagerTest.getAllFromSQLTable("sql/unigene.db", "unigene");
        SQLInputStream is = new SQLInputStream(rs, false, 1, 3);

        FeatureCodec featCodec = new ExampleBinaryCodec();
        PositionalBufferedStream pbs = new PositionalBufferedStream(is);
        int count = 0;
        while (!rs.isAfterLast()) {
            Feature feat = featCodec.decode(pbs);
            checkFeatureIntegrity(feat, "chr2");
            count++;
        }
        assertEquals(72, count);

        DBManager.closeAll(rs);
    }

    @Test
    public void basicTestReadString() throws Exception {
        ResultSet rs = DBManagerTest.getAllFromSQLTable("sql/unigene.db", "unigene");
        InputStream is = new SQLInputStream(rs, true);
        BufferedReader in = new BufferedReader(new InputStreamReader(is));
        int totalLines = 0;
        String val = "";
        AsciiFeatureCodec featureCodec = new IGVBEDCodec(genome);
        while ((val = in.readLine()) != null) {
            Feature feat = featureCodec.decode(val);
            checkFeatureIntegrity(feat, "chr2");
            totalLines += 1;
        }
        assertEquals(72, totalLines);

        DBManager.closeAll(rs);
    }

    private void checkFeatureIntegrity(Feature feat, String expChr) throws Exception {
        assertEquals(expChr, feat.getChr());
        assertTrue(feat.getEnd() > feat.getStart());
        assertTrue(feat.getEnd() > 0);
        assertTrue(feat.getStart() > 0);
    }
}
