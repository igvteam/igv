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

package org.broad.igv.tools.sort;

import com.google.common.base.Function;
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

import static org.junit.Assert.*;


/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Sep 18, 2009
 * Time: 6:19:28 PM
 */
public class SorterTest extends AbstractHeadlessTest {

    @Test
    public void testSortBed() throws Exception {
        testSort(TestUtils.DATA_DIR + "bed/Unigene.unsorted.bed", 0, 1, 10, 71, 0);
    }

    @Test
    public void testSortBed1() throws Exception {
        testSort(TestUtils.DATA_DIR + "bed/Unigene.unsorted1.bed", 0, 1, 10, 71, 0);
    }

    @Test
    public void testSortBed2() throws Exception {
        testSort(TestUtils.DATA_DIR + "bed/GSM1004654_10k.bed", 0, 1, 50, 10000, 0);
    }

    @Test
    public void testSortVCF() throws Exception {
        testSort(TestUtils.DATA_DIR + "vcf/SRP32_v4.0.vcf", 0, 1);
    }

    @Test
    public void testSortGFF() throws Exception {
        testSort(TestUtils.DATA_DIR + "gff/aliased.unsorted.gff", 0, 3);
    }

    @Test
    public void testSortCN() throws Exception{
        String path = TestUtils.DATA_DIR + "cn/1klines.cn";
        testSort(path, 1, 2, 10, 1000, 1);
    }

    @Test
    public void testSortGWAS() throws Exception{
        String path = TestUtils.DATA_DIR + "gwas/random.gwas";
        testSort(path, 0, 1, 10, 100, 1);
    }

    public void testSort(String infile, int chrCol, int startCol) throws IOException {
        testSort(infile, chrCol, startCol, 10, null, 0);
    }

    public void testSort(String infile, int chrCol, int startCol, int maxRecords, Integer expectedLines, int skipTopLines) throws IOException {

        File ifile = new File(infile);
        File ofile = new File(infile + ".sorted");
        ofile.deleteOnExit();

        Sorter sorter = Sorter.getSorter(ifile, ofile);
        sorter.setMaxRecords(maxRecords);
        sorter.run();

        int outLines = checkFileSorted(ofile, chrCol, startCol, skipTopLines);

        if(expectedLines != null){
            assertEquals((int) expectedLines, outLines);
        }
    }

    public static int checkFileSorted(File ofile, int chrCol, int startCol, int skipTopLines) {
        BufferedReader reader = null;
        int numlines = 0;
        String nextLine = "";

        try {
            reader = new BufferedReader(new FileReader(ofile));
            String lastChr = "";
            int lastStart = 0;
            Set<String> chromosomes = new HashSet();
            for(int ii=0; ii < skipTopLines; ii++){
                reader.readLine();
            }
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
                    String msg = String.format("Chromosome %s out of order in line: %s", chr, nextLine);
                    assertFalse(msg, chromosomes.contains(chr));
                    chromosomes.add(chr);
                }
                numlines++;
                lastChr = chr;
                lastStart = start;
            }
        } catch (Exception e) {
            e.printStackTrace();
            String msg = "Exception during checking on line: " + nextLine;
            throw new AssertionError(msg);
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

        Function<Sorter, Void> func = new Function<Sorter, Void>() {
            @Override
            public Void apply(Sorter input) {
                try {
                    input.run();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return null;
            }
        };

        long[] times = TestUtils.timeMethod(supplier, func, nTrials);
    }

}
