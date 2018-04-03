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
