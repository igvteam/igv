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

package org.broad.igv.tools.sort;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;


/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Sep 18, 2009
 * Time: 6:19:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class SorterTest {


    @Test
    public void testSortBed() throws Exception {
        testSort(TestUtils.DATA_DIR + "/bed/Unigene.unsorted.bed");
    }

    @Test
    public void testSortVCF() throws Exception {
        testSort(TestUtils.DATA_DIR + "/vcf/SRP32_v4.0.vcf");
    }

    public void testSort(String infile) throws IOException {

        File ifile = new File(infile);
        File ofile = new File(infile + ".sorted");
        ofile.deleteOnExit();

        Sorter sorter = Sorter.getSorter(ifile, ofile);
        sorter.setMaxRecords(10);  // <= force text of serialization
        sorter.run();

        checkBedSorted(ofile);
    }

    public static int checkBedSorted(File ofile) {
        BufferedReader reader = null;
        int numlines = 0;
        try {
            reader = new BufferedReader(new FileReader(ofile));
            String nextLine = "";
            String lastChr = "";
            int lastStart = 0;
            Set<String> chromosomes = new HashSet();
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("track") || nextLine.startsWith(("#"))) {
                    continue;
                }
                String[] tokens = nextLine.split("\t");
                String chr = tokens[0];

                int start = Integer.parseInt(tokens[1]);
                if (chr.equals(lastChr)) {
                    assertTrue(start >= lastStart);
                } else {
                    assertFalse(chromosomes.contains(chr));
                    chromosomes.add(chr);
                }
                numlines++;
                lastChr = chr;
                lastStart = start;
            }
        } catch (Exception e) {
            throw new AssertionError("Exception during checking: " + e);
        } finally {
            try {
                reader.close();
            } catch (Exception e) {

            }
        }
        return numlines;
    }
}
