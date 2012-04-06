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
package org.broad.igv.feature;

import org.broad.igv.util.TestUtils;
import org.junit.*;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.fail;

/**
 * @author jrobinso
 */
public class ProbeToGeneMapTest {

    public ProbeToGeneMapTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of affymetrix critera.
     */
    @Test
    public void testAffy() throws FileNotFoundException, IOException {
        File f = new File(TestUtils.LARGE_DATA_DIR + "affy_probe_gene_mapping.txt");
        BufferedReader br = new BufferedReader(new FileReader(f));
        String nextLine = br.readLine();
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            final String probe = tokens[0];
            if (!isAffy(probe)) {
                fail("Probe: " + probe + " failed Affy test criteria.");
            } else {
                // This is an affy probe by our test.  See if it matches other platforms.
                if (isAgilent(probe)) {
                    fail("Agilent conflict: " + probe);
                }
                if (isIllumina(probe)) {
                    fail("Illumina conflict: " + probe);
                }
            }
        }
        br.close();
    }

    @Test
    public void testAgilent() throws FileNotFoundException, IOException {
        File f = new File(TestUtils.LARGE_DATA_DIR + "agilent_probe_gene_mapping.txt");
        BufferedReader br = new BufferedReader(new FileReader(f));
        String nextLine = br.readLine();

        // See if there is any overlap with affy test      
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            final String probe = tokens[0];
            if (!isAgilent(probe)) {
                fail("Probe: " + probe + " failed Agilent test criteria.");
            } else {
                // This is an agilent probe by our test.  See if it matches other platforms.
                if (isAffy(probe)) {
                    fail("Affy conflict: " + probe);
                }
                if (isIllumina(probe)) {
                    fail("Illumina conflict: " + probe);
                }
            }
        }
        br.close();
    }

    @Test
    public void testIllumina() throws FileNotFoundException, IOException {
        File f = new File(TestUtils.LARGE_DATA_DIR + "illumina_probe_gene_mapping.txt");
        BufferedReader br = new BufferedReader(new FileReader(f));
        String nextLine = br.readLine();

        // See if there is any overlap with affy test   
        List<String> noMatchList = new ArrayList();
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            final String probe = tokens[0];
            if (!isIllumina(probe)) {
                noMatchList.add(probe);
                //fail("Probe: " + probe + " failed Illumina test criteria.");
            } else {
                // This is an agilent probe by our test.  See if it matches other platforms.
                if (isAffy(probe)) {
                    fail("Affy conflict: " + probe);
                }
                if (isAgilent(probe)) {
                    fail("Agilent conflict: " + probe);
                }
            }
        }
        System.out.println("Non matching illumina probes:");

        for (String probe : noMatchList) {
            System.out.println(probe);
        }

        br.close();
    }

    private boolean isIllumina(String probe) {
        return probe.startsWith("ILMN_") || probe.startsWith("GI_") || probe.startsWith("NM_") ||
                probe.startsWith("XM_");
    }

    private boolean isAffy(String probe) {
        return probe.endsWith("_at") || probe.endsWith("_st");
    }

    private boolean isAgilent(String probe) {
        return probe.startsWith("A_");
    }

    /**
     * Test of getGenesForProbe method, of class ProbeToLocusMap.
     */
    @Test
    public void getGenesForProbe() {
        // TODO -- implementation
    }
}
