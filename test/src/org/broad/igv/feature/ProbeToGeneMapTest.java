/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.junit.*;
import static org.junit.Assert.fail;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

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
        File f = new File("../server/web/resources/probes/affy_probe_gene_mapping.txt");
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
        File f = new File("../server/web/resources/probes/agilent_probe_gene_mapping.txt");
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
        File f = new File("../server/web/resources/probes/illumina_probe_gene_mapping.txt");
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
