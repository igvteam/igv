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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.reader.BAMWebserviceReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
@Ignore
public class BAMRemoteQueryReaderTest extends AbstractHeadlessTest {

    public BAMRemoteQueryReaderTest() {
    }

    static PreferenceManager preferenceManager;

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        Globals.READ_TIMEOUT = 5 * 30 * 1000;
        Globals.CONNECT_TIMEOUT = 5 * 30 * 1000;
        preferenceManager = PreferenceManager.getInstance();
    }

    /**
     * Test of query method, of class BAMRemoteQueryReader.
     */
    @Test
    public void testQuery() throws IOException {
        /*
        http://www.broadinstitute.org/webservices/igv/bam?method=query&samFile=/broad/1KG/DCC_merged/freeze5/NA12878.pilot2.SLX.bam&chr=16&start=50542554&end=50542722
        String path = "/Users/jrobinso/IGV/Alignments/303KY.8.paired.bam";
        String serverURL = "http://localhost:8080/webservices/igv";
        String chr = "chr7";
        int start = 2244149;
        int end = 2250144;
         * */

        String path = "/1KG/prod/data/HG00099/alignment/HG00099.chrom20.SOLID.bfast.GBR.low_coverage.20101123.bam";
        String serverURL = "http://www.broadinstitute.org/igvdata/";
        String fullpath = serverURL + path;
        String chr = "chr1";
        int start = 50542554; //557000;  //
        int end = 50542722; //558000; //
        boolean contained = false;


        //ResourceLocator rl = new ResourceLocator(fullpath);
        //checkReader(rl, chr, start, end, contained);

        ResourceLocator rlremote = new ResourceLocator(serverURL, path);
        checkReader(rlremote, chr, start, end, contained);
    }

    @Test
    public void testQuery2(){
        String path = "/1KG/DCC_merged/freeze4/NA12878.ceu.daughter.bam";
        String serverURL = "http://www.broadinstitute.org/igvdata/";

        String chr = "1";
        int start = 713700;
        int end = 714100;

        ResourceLocator locator = new ResourceLocator(serverURL, path);
        checkReader(locator, chr, start, end, true);

    }

    private void checkReader(ResourceLocator locator, String chr, int start, int end, boolean contained){
        Alignment al = null;
        BAMWebserviceReader instance = new BAMWebserviceReader(locator);
        CloseableIterator<Alignment> result = instance.query(chr, start, end, contained);
        //long t0 = System.currentTimeMillis();
        int count = TestUtils.assertFeatureIteratorSorted(result);

        assertTrue("No data received", count > 0);

    }
}