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
package org.broad.igv.util;

import java.io.*;
import java.util.HashSet;

/**
 * @author jrobinso
 */
public class IlluminaUtils {

    public static final String LOCUS_DELIMITER = "///";

    public static void main(String[] args) throws IOException {
        File annotFile = new File(
                "/Users/jrobinso/IGV/ProbeAnnotations/illumina");
        File outputFile = new File(
                "/Users/jrobinso/IGV/ProbeAnnotations/illumina_probe_gene_mapping.txt");
        createProbeMappingFile(annotFile, outputFile);
    }

    static void createProbeMappingFile(File inputDir, File outputFile) throws IOException {
        BufferedReader reader = null;
        PrintWriter pw = null;
        HashSet<String> processedProbes = new HashSet(100000);
        try {
            int dupCount = 0;
            pw = new PrintWriter(outputFile);
            for (File inputFile : inputDir.listFiles()) {
                System.out.println(inputFile.getName());
                reader = new BufferedReader(new FileReader(inputFile));
                String nextLine;
                while ((nextLine = reader.readLine()) != null) {
                    String[] tokens = nextLine.split("\t");
                    boolean isMIArray = inputFile.getName().endsWith("MAP.txt");
                    int probeColumn = isMIArray ? 3 : 13;
                    int chrColumn = isMIArray ? 19 : 18;
                    int coordColumn = 20;
                    int probeSeqColumn = 4;

                    if (tokens.length > coordColumn) {

                        int probeLength = tokens[probeSeqColumn].trim().length();

                        String probe = tokens[probeColumn].trim();
                        if (processedProbes.contains(probe)) {
                            dupCount++;
                        } else {
                            processedProbes.add(probe);

                            String chr = "chr" + tokens[chrColumn].trim();
                            String[] coords = tokens[coordColumn].split(":");

                            pw.print(probe + "\t");
                            for (int i = 0; i < coords.length; i++) {
                                pw.print(chr + ":" + coords[i]);
                                if (isMIArray) {
                                    pw.print("-" + (Integer.parseInt(coords[i]) + probeLength));
                                }
                                if (i < coords.length - 1) {
                                    pw.print(LOCUS_DELIMITER);
                                }
                            }
                            pw.println();
                        }
                    }
                }
                reader.close();
                reader = null;
            }
            System.out.println("# overlapping probes = " + dupCount);
        } finally {
            if (reader != null) {
                reader.close();
            }
            pw.close();
        }

    }
}
