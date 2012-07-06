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
package org.broad.igv.tools;

import org.broad.igv.util.HttpUtils;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.zip.GZIPInputStream;

/**
 * @author jrobinso
 */
public class TestFileGenerator {

    public static void main(String[] args) {

        generateTestFile("/Users/jrobinso/IGV/test_25thousand.gct", true, 25000, 100);
    }

    enum FileType {

        CN, IGV, GCT, WIG
    }

    ;

    public static void generateTestFile(
            String outputFile,
            boolean sorted,
            int nRows,
            int nSamples) {

        int rowsPerChr = Math.max(1, nRows / chromosome.length);

        PrintWriter pw = null;

        FileType type = null;
        if (outputFile.endsWith(".cn")) {
            type = FileType.CN;
        } else if (outputFile.endsWith(".igv")) {
            type = FileType.IGV;
        } else if (outputFile.endsWith(".gct")) {
            type = FileType.GCT;

        } else if (outputFile.endsWith(".wig")) {
            type = FileType.WIG;

        } else {
            System.out.println("Unsupported file type: " + outputFile);
        }

        try {
            pw = new PrintWriter(new FileWriter(outputFile));
            writeHeader(pw, nSamples, type);

            Random rand = new Random();
            List<String> probes = type == FileType.GCT ? getExpressionProbes() : new ArrayList();
            int n = 0;

            for (int i = 0; i < chromosome.length; i++) {
                String chr = chromosome[i];
                int len = chromSize[i];
                double delta = ((double) len) / rowsPerChr;
                double center = -3;
                int lastStart = 0;

                if (type == FileType.WIG) {
                    int step = (int) delta;
                    pw.println("fixedStep chrom=" + chr + " start=" + 1 + " step=" + step + " span=" + (step - 1));
                }

                for (int r = 0; r <= rowsPerChr; r++) {

                    int start = getStart(len, lastStart, r, delta, rand, sorted);
                    lastStart = start + 1;

                    center = 3 * ((double) start / len);

                    if (type == FileType.WIG) {
                        pw.println(center + rand.nextGaussian());

                    } else {
                        switch (type) {
                            case CN:
                                pw.print("Snp_" + r + "\t" + chr + "\t" + start);
                                break;
                            case IGV:
                                int end = Math.max(start + 1, (int) (start + Math.random() * 1000));
                                pw.print(chr + "\t" + start + "\t" + end + "\tProbe_" + r);
                                break;
                            case GCT:
                                String probe = probes.get(n);
                                pw.print(probe + "\tdesc_" + n);
                        }

                        for (int s = 0; s < nSamples; s++) {
                            double v = center + rand.nextGaussian();
                            pw.print("\t" + v);
                        }

                        pw.println();
                        n++;
                        if (n >= probes.size()) {
                            n = 0;
                        }
                    }
                }
            }
        } catch (IOException e) {
            System.out.println("Error: " + e.getMessage());
        } finally

        {
            pw.close();
        }

    }

    public static int getChromosomeIndex(int i, boolean sorted) {
        if (sorted) {
            return i;
        } else {
            return (int) (Math.random() * 24);
        }
    }

    public static int getStart(int len, int lastStart, int r, double delta, Random rand, boolean sorted) {
        if (sorted) {
            return Math.min(len, Math.max(lastStart, (int) (r * (delta - 1) - rand.nextGaussian() * delta)));
        } else {
            return (int) (Math.random() * len);
        }

    }

    public static void writeHeader(PrintWriter pw, int nSamples, FileType type) {

        if (type == FileType.CN) {
            pw.print("Snp\tChr\tPos");
            for (int s = 0; s < nSamples; s++) {
                pw.print("\tSample_" + s);
            }
            pw.println();
        } else if (type == FileType.IGV) {
            pw.print("Chr\tStart\tEnd\tProbe");
            for (int s = 0; s < nSamples; s++) {
                pw.print("\tSample_" + s);
            }
            pw.println();
        } else if (type == FileType.GCT) {
            pw.println("ignored");
            pw.println("ignored");
            pw.print("Probe\tdesc");
            for (int s = 0; s < nSamples; s++) {
                pw.print("\tSample_" + s);
            }
            pw.println();
        }
    }

    static List<String> getExpressionProbes() {
        String urlString = "http://www.broadinstitute.org/igv/resources/probes/affy/affy_human_mappings.txt.gz";
        AsciiLineReader bufReader = null;
        List<String> probes = new ArrayList(60000);
        try {
            InputStream is = null;

            if (HttpUtils.isRemoteURL(urlString)) {
                URL url = new URL(urlString);
                is = HttpUtils.getInstance().openConnectionStream(url);
            } else {
                is = new FileInputStream(urlString);
            }
            if (urlString.endsWith("gz")) {
                is = new GZIPInputStream(is);
            }
            bufReader = new AsciiLineReader(is);

            String nextLine = "";
            while ((nextLine = bufReader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                probes.add(tokens[0]);
            }

            return probes;

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        } finally {
            if (bufReader != null) {
                bufReader.close();
            }
        }

    }

    static String[] chromosome = {
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            "chrX",
            "chrY",
            "chrM"};
    static int[] chromSize = {
            247249719,
            242951149,
            199501827,
            191273063,
            180857866,
            170899992,
            158821424,
            146274826,
            140273252,
            135374737,
            134452384,
            132349534,
            114142980,
            106368585,
            100338915,
            88827254,
            78774742,
            76117153,
            63811651,
            62435964,
            46944323,
            49691432,
            154913754,
            57772954,
            16571};
}
