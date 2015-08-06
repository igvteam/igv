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
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.tribble.UCSCGeneTableCodec;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import static junit.framework.Assert.*;

/**
 * Test our profiles contacting UCSC database
 * User: jacob
 * Date: 2012/05/29
 */
public class UCSC_SQL_Test extends AbstractHeadlessTest {

    private static final String UCSC_HOST = "genome-mysql.cse.ucsc.edu";
    private String profilePath = TestUtils.DATA_DIR + "sql/UCSC_profiles.dbxml";

    private static void checkConnectUCSC() throws Exception {
        boolean succeeded = false;
        try {
            String testPath = DBManager.createConnectionURL("mysql", UCSC_HOST, "hg18", null);
            ResourceLocator locator = new ResourceLocator(testPath);
            locator.setUsername("genome");
            DBManager.getConnection(locator);
            succeeded = true;
        } catch (Exception e) {
            e.printStackTrace();
        }

        Assume.assumeTrue(succeeded);
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        checkConnectUCSC();
    }

    @Test
    public void testLoadUCSC() throws Exception {
        AsciiFeatureCodec codec = new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.UCSCGENE, genome);

        String host = UCSC_HOST;

        String path = "hg18";
        String port = null;

        String url = DBManager.createConnectionURL("mysql", host, path, port);
        ResourceLocator locator = new ResourceLocator(url);
        locator.setUsername("genome");

        String tableName = "knownGene";
        int strt = 100000;
        int end = 400000;

        DBProfile.DBTable table = new DBProfile.DBTable(locator, tableName, "n/a", null, SQLCodecSource.UCSC_CHROMO_COL, SQLCodecSource.UCSC_START_COL, SQLCodecSource.UCSC_END_COL, 1, Integer.MAX_VALUE, null, null, null);
        SQLCodecSource reader = new SQLCodecSource(table, codec);
        Iterator<Feature> SQLFeatures = reader.getFeatures("chr1", strt, end);

        int count = 0;
        while (SQLFeatures.hasNext()) {
            SQLFeatures.next();
            count++;
        }
        assertEquals(12, count);

        /**
         * We should be getting unique names only
         */
        List<String> names = reader.getSequenceNames();
        assertTrue(names.size() > 0);

        Set<String> set = new HashSet<String>(names);
        assertEquals(names.size(), set.size());

    }

    @Test
    public void testLoadUCSCFromProfileGene() throws Exception {
        tstLoadFromProfile(profilePath, "knownGene");
    }

    @Test
    public void testLoadUCSCFromProfileBED() throws Exception {
        tstLoadFromProfile(profilePath, "affyExonTissues");
    }

    @Test
    public void testLoadUCSCFromProfilePSL() throws Exception {
        tstLoadFromProfile(profilePath, "all_mrna");
    }

    @Test
    public void testLoadUCSCFromProfileSNP() throws Exception {
        SQLCodecSource source = tstLoadFromProfile(profilePath, "snp126");
        Iterator<Feature> feats = source.getFeatures("chr2", 10000, 100000);
        int count = 0;
        while (feats.hasNext()) {
            BasicFeature f = (BasicFeature) feats.next();
            assertEquals(0.0f, f.getScore());
            assertFalse(f.hasExons());
            assertNotSame(Strand.NONE, f.getStrand());
            count++;
        }

        assertTrue(count > 0);
    }

    public SQLCodecSource tstLoadFromProfile(String profilePath, String tableName) throws Exception {
        SQLCodecSource source = SQLCodecSource.getFromProfile(profilePath, tableName);
        if (source == null) {
            throw new RuntimeException("Table " + tableName + " not found in profile " + profilePath);
        }
        int start = 1;
        int end = 100000;
        Iterator<Feature> feats = source.getFeatures("chr1", start, end);
        int count = 0;

        while (feats.hasNext()) {
            Feature f = feats.next();
            assertTrue(f.getStart() >= start);
            assertTrue(f.getStart() < end);
            count++;
        }

        assertTrue("No data retrieved", count > 0);

        return source;

    }

    @Test
    public void testQueryWithBins() throws Exception {
        tstQueryWithBins(profilePath, "all_mrna", "chr3", 500, 500000);
    }

    @Ignore("Fails sometimes and we're not sure why, but functionality isn't supported yet")
    @Test
    public void testQueryWithBins_big() throws Exception {
        tstQueryWithBins(profilePath, "all_mrna", "chr1", 0, (int) 247e4);
    }

    public void tstQueryWithBins(String profilePath, String tableName, String chr, int start, int end) throws Exception {
        int featWindowSize = (end - start) * 2;
        SQLCodecSource binSource = SQLCodecSource.getFromProfile(profilePath, tableName);
        binSource.setFeatureWindowSize(featWindowSize);
        assertNotNull(binSource.binColName);

        SQLCodecSource noBinSource = SQLCodecSource.getFromProfile(profilePath, tableName);
        noBinSource.setFeatureWindowSize(featWindowSize);
        binSource.binColName = null;


        Iterator<Feature> binFeats = binSource.getFeatures(chr, start, end);
        Iterator<Feature> noBinFeats = noBinSource.getFeatures(chr, start, end);
        assertNotNull(binFeats);
        assertNotNull(noBinFeats);

        int count = 0;
        while (binFeats.hasNext()) {
            assertTrue(noBinFeats.hasNext());

            Feature bf = binFeats.next();
            Feature nbf = noBinFeats.next();

            assertEquals(chr, bf.getChr());
            assertTrue(bf.getStart() >= start);
            assertTrue(bf.getStart() < end);

            assertEquals(nbf.getChr(), bf.getChr());
            assertEquals(nbf.getStart(), bf.getStart());
            assertEquals(nbf.getEnd(), bf.getEnd());

            count++;
        }

        assertFalse(noBinFeats.hasNext());


    }
}
