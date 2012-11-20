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

import com.google.common.base.Predicate;
import com.google.common.base.Supplier;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.genome.ChromosomeNameComparator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Comparator;
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
public class SorterTest extends AbstractHeadlessTest {

    @Test
    public void testSortBed() throws Exception {
        testSort(TestUtils.DATA_DIR + "bed/Unigene.unsorted.bed", 0, 1);
    }

    @Test
    public void testSortBed1() throws Exception {
        testSort(TestUtils.DATA_DIR + "bed/Unigene.unsorted1.bed", 0, 1);
    }

    @Test
    public void testSortBed2() throws Exception {
        testSort(TestUtils.DATA_DIR + "bed/GSM1004654_10k.bed", 0, 1, 50);
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
        testSort(infile, chrCol, startCol, 10);
    }

    public void testSort(String infile, int chrCol, int startCol, int maxRecords) throws IOException {

        File ifile = new File(infile);
        File ofile = new File(infile + ".sorted");
        ofile.deleteOnExit();

        Sorter sorter = Sorter.getSorter(ifile, ofile);
        sorter.setMaxRecords(maxRecords);
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


    //@Test
    public void testCurrentCompSpeed() throws IOException {
        tstComparatorSpeed(Sorter.getDefaultComparator());
    }

    //@Test
//    public void testTestCompSpeed() throws IOException{
//        tstComparatorSpeed(getTestComparator());
//    }

    public void tstComparatorSpeed(final Comparator<SortableRecord> testComp) throws IOException {
        final File inputFile = new File(TestUtils.DATA_DIR + "bed", "GSM1004654_10k.bed");
        final File outputFile = new File(TestUtils.TMP_OUTPUT_DIR, "GSM1004654_10k.sorted.bed");

        int nTrials = 100;

        Supplier<Sorter> supplier = new Supplier<Sorter>() {
            @Override
            public Sorter get() {
                ChromosomeNameComparator.get().resetCache();
                Sorter sorter = new BedSorter(inputFile, outputFile);
                sorter.setComparator(testComp);
                return sorter;
            }
        };

        Predicate<Sorter> predicate = new Predicate<Sorter>() {
            @Override
            public boolean apply(Sorter input) {
                try {
                    input.run();
                } catch (IOException e) {
                    return false;
                }
                return true;
            }
        };

        long[] times = TestUtils.timeMethod(supplier, predicate, nTrials);
    }

}
