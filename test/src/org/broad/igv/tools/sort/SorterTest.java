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
 */
public class SorterTest {


    @Test
    public void testSortBed() throws Exception {
        testSort(TestUtils.DATA_DIR + "bed/Unigene.unsorted.bed", 0, 1);
    }

    @Test
    public void testSortVCF() throws Exception {
        testSort(TestUtils.DATA_DIR + "vcf/SRP32_v4.0.vcf", 0, 1);
    }

    @Test
    public void testSortGFF() throws Exception {
        testSort(TestUtils.DATA_DIR + "gff/aliased.unsorted.gff", 0, 3);
    }

    public void testSort(String infile, int chrCol, int startCol) throws IOException {

        File ifile = new File(infile);
        File ofile = new File(infile + ".sorted");
        ofile.deleteOnExit();

        Sorter sorter = Sorter.getSorter(ifile, ofile);
        sorter.setMaxRecords(10);  // <= force text of serialization
        sorter.run();

        checkFileSorted(ofile, chrCol, startCol);
    }

    public static int checkFileSorted(File ofile, int chrCol, int startCol) {
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
                String chr = tokens[chrCol];

                int start = Integer.parseInt(tokens[startCol]);
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
