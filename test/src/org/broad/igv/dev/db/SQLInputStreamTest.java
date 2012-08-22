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

package org.broad.igv.dev.db;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.example.ExampleBinaryCodec;
import org.broad.tribble.readers.PositionalBufferedStream;
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
    }

    private void checkFeatureIntegrity(Feature feat, String expChr) throws Exception {
        assertEquals(expChr, feat.getChr());
        assertTrue(feat.getEnd() > feat.getStart());
        assertTrue(feat.getEnd() > 0);
        assertTrue(feat.getStart() > 0);
    }
}
