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
import org.broad.igv.feature.tribble.UCSCGeneTableCodec;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.junit.Test;

import java.io.File;
import java.util.Iterator;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012/05/29
 */
public class SQLCodecReaderTest extends AbstractHeadlessTest {

    @Test
    public void testLoadBED() throws Exception {

        FeatureCodec codec = new IGVBEDCodec();
        SQLCodecReader reader = new SQLCodecReader(null, codec);

        ResourceLocator locator = new ResourceLocator(null);
        String host = (new File(TestUtils.DATA_DIR)).getAbsolutePath();
        String path = "sql/unigene.db";
        locator.setUrl(DBManager.createConnectionURL(host, path, null,"sqlite"));
        locator.setDescription("SELECT * FROM unigene ORDER BY chrom, chromStart");
        Iterable<Feature> SQLFeatures = reader.load(locator, null);

        String bedFile = host + "/bed/unigene.sample.bed";
        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, codec, false);
        Iterator<Feature> fileFeatures = bfr.iterator();

        int count=0;
        for(Feature f: SQLFeatures){
            Feature fileFeature = fileFeatures.next();
            assertEquals(fileFeature.getChr(), f.getChr());
            assertEquals(fileFeature.getStart(), f.getStart());
            assertEquals(fileFeature.getEnd(), f.getEnd());
            count++;
        }

        assertEquals(72, count);
    }

    @Test
    public void testLoadUCSC() throws Exception{
        FeatureCodec codec = new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.UCSCGENE, genome);
        SQLCodecReader reader = new SQLCodecReader(null, codec);

        ResourceLocator locator = new ResourceLocator(null);

        String host = "genome-mysql.cse.ucsc.edu";

        String path = "hg19";
        String port = null;

        DBManager.username = "genome";

        locator.setUrl(DBManager.createConnectionURL(host, path, port, "mysql"));
        String table = "knownGene";
        int strt = 100000;
        int end = 400000;
        String query = String.format("SELECT * FROM %s WHERE chrom = 'chr1' AND txStart >= %d AND txStart < %d ORDER BY txStart;", table, strt, end );

        locator.setDescription(query);
        Iterable<Feature> SQLFeatures = reader.load(locator, null);

        int count=0;
        for(Feature f: SQLFeatures){
            count++;
        }

        assertEquals(9, count);





    }
}
