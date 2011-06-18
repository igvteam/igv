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
package org.broad.igv.feature;

//~--- JDK imports ------------------------------------------------------------

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * chromosome
 * coding_exon
 * contig
 * EST
 * EST_match
 * exon
 * expressed_sequence_match
 * gap
 * gene
 * ##gff-version 3
 * match
 * microarray_oligo
 * mRNA
 * processed_transcript
 * region
 *
 * @author jrobinso
 */
public class FeatureFileUtils {


    public static void main(String[] args) throws IOException {
        //String gffFile = "/Users/jrobinso/celegans/current.gff3";
        //String outputDir = "/Users/jrobinso/celegans";
        //splitGFFFileByType(gffFile, outputDir);
        covertProbeMapToBedFile("/Users/jrobinso/IGV/TestData/expression/HG-U133_PLUS_2.mapping.txt",
                "/Users/jrobinso/IGV/TestData/expression/HG-U133_PLUS_2.mapping.bed");
    }

    static void covertProbeMapToBedFile(String probeMapFile, String bedFile) throws
            FileNotFoundException, IOException {


        BufferedReader br = new BufferedReader(new FileReader(probeMapFile));
        PrintWriter pw = new PrintWriter(new FileWriter(bedFile));

        String nextLine;
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            Locus locus = new Locus(tokens[1].trim());
            pw.println(
                    locus.getChr() + "\t" + locus.getStart() + "\t" + locus.getEnd() + "\t" + tokens[0].trim());
        }

        br.close();
        pw.close();

    }

    static void splitEmblFileByType(String emblFile, String outputDirectory) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(emblFile));
        String nextLine;
        Set<String> codes = new HashSet();

        while ((nextLine = br.readLine()) != null) {
            if (nextLine.startsWith("FT") && (nextLine.length() > 19)) {
                String code = nextLine.substring(5, 19).trim();
                if (code.length() > 0) {
                    codes.add(code);
                }
            }
        }
        br.close();

        Map<String, PrintWriter> writers = new HashMap();
        for (String code : codes) {
            writers.put(code,
                    new PrintWriter(new FileWriter(new File(outputDirectory, code + ".embl"))));
        }

        br = new BufferedReader(new FileReader(emblFile));
        PrintWriter currentWriter = null;
        while ((nextLine = br.readLine()) != null) {
            if (nextLine.startsWith("ID")) {
                for (PrintWriter pw : writers.values()) {
                    pw.println(nextLine);
                }
            } else if (nextLine.startsWith("FT")) {
                String code = nextLine.substring(5, 19).trim();
                if (code.length() > 0) {
                    currentWriter = writers.get(code);
                }
                if (currentWriter != null) {
                    currentWriter.println(nextLine);
                }
            } else {
                currentWriter = null;
            }
        }

        br.close();
        for (PrintWriter pw : writers.values()) {
            pw.close();
        }
    }
}
