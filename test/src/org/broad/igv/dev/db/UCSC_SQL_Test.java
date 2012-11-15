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
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.tribble.UCSCGeneTableCodec;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.junit.Assume;
import org.junit.BeforeClass;
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

        DBTable table = new DBTable(locator, tableName, "n/a", null, SQLCodecSource.UCSC_CHROMO_COL, SQLCodecSource.UCSC_START_COL, SQLCodecSource.UCSC_END_COL, 1, Integer.MAX_VALUE, null, null, null);
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

//        for(int att=0; att < 10; att++){
//            source.queryStatement.setString(1, "fake");
//            source.queryStatement.setInt(3, 93);
//            Iterator<Feature> feats2 = source.getFeatures("chr1", start, end);
//            int tCount = 0;
//            while(feats2.hasNext()){
//                tCount++;
//            }
//            assertEquals(count, tCount);
//        }

        assertNotNull(source.queryStatement);

        return source;

    }

    @Test
    public void testQueryWithBins() throws Exception {
        tstQueryWithBins(profilePath, "all_mrna", "chr3", 500, 500000);
    }

    @Test
    public void testQueryWithBins_big() throws Exception {
        tstQueryWithBins(profilePath, "all_mrna", "chr1", 0, (int) 247e4);
    }

    public void tstQueryWithBins(String profilePath, String tableName, String chr, int start, int end) throws Exception {
        SQLCodecSource binSource = SQLCodecSource.getFromProfile(profilePath, tableName);
        assertNotNull(binSource.binColName);

        SQLCodecSource noBinSource = SQLCodecSource.getFromProfile(profilePath, tableName);
        binSource.binColName = null;


        Iterator<Feature> binFeats = binSource.getFeatures(chr, start, end);
        Iterator<Feature> noBinFeats = noBinSource.getFeatures(chr, start, end);

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
