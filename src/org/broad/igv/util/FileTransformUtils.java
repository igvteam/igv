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

package org.broad.igv.util;

import org.broad.igv.Globals;

import java.io.*;
import java.util.*;

/**
 * Mostly one-off static methods to transform or maniuplate feature file formats.
 */
public class FileTransformUtils {

    static Set<String> types = new HashSet(Arrays.asList("SINE", "LINE", "LTR", "DNA", "Simple_repeat",
            "Low_complexity", "Satellite", "RNA", "Other", "Unknown", "Uncategorized"));

    public static void main(String[] args) throws IOException {
        String ifile = "/Users/jrobinso/projects/IGV/humanv3_hg18Pos.csv";
        String output = "/Users/jrobinso/projects/IGV/humanv3Pos_hg18.bed";
        probeToBed(ifile, output, true);
    }

    /**
     * Utility to convert a file with the following format to "bed"
     * ILMN_2087817	chrY:9529429:9529478:+					        // "1" based?
     * ILMN_2204360	chrY:9598604:9598615:-	chrY:9590985:9591022:-
     *
     * @param iFile
     * @param oFile
     * @param includeMultiMappings if false probes with multiple location mappings are filtered out
     */
    public static void probeToBed(String iFile, String oFile, boolean includeMultiMappings) throws IOException {

        BufferedReader br = null;
        PrintWriter pw = null;

        try {
            br = ParsingUtils.openBufferedReader(iFile);
            pw = new PrintWriter(new FileWriter(oFile));

            String nextLine;
            br.readLine();  // eat header
            while ((nextLine = br.readLine()) != null) {

                String[] tokens = Globals.commaPattern.split(nextLine);
                String probe = tokens[0];

                if(tokens.length < 1 || (tokens.length > 2 && !includeMultiMappings)) continue;

                for (int i = 1; i < tokens.length; i++) {

                    String loc = tokens[i];
                    String[] locParts = Globals.colonPattern.split(loc);
                    String chr = locParts[0];
                    int start = Integer.parseInt(locParts[1]) - 1;
                    int end = Integer.parseInt(locParts[2]);
                    String strand = locParts[3];
                    pw.println(chr + "\t" + start + "\t" + end + "\t" + probe + "\t1000\t" + strand);
                }
            }
        } finally {
            if(br != null) br.close();
            if(pw != null) pw.close();
        }


    }

    public static void splitRepeatMasker(String iFile, String outputDirectory, String prefix) throws IOException {


        BufferedReader br = null;
        Map<String, PrintWriter> pws = new HashMap();

        try {
            br = new BufferedReader(new FileReader(iFile));
            File dir = new File(outputDirectory);
            for (String type : types) {
                File f = new File(dir, prefix + type + ".bed");
                pws.put(type, new PrintWriter(new BufferedWriter(new FileWriter(f))));
            }

            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("#")) continue;
                String[] tokens = nextLine.split("\t");
                String type = getType(tokens[5]);

                // Rerrange columns for legal bed
                pws.get(type).println(tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] + "\t" +
                        tokens[4] + "\t" + tokens[3]);
            }


        } finally {
            if (br != null) {
                br.close();
            }
            for (PrintWriter pw : pws.values()) {
                pw.close();
            }
        }

    }


    public static String getType(String s) {

        s = s.replace("?", "");

        if (s.contains("RNA")) {
            return "RNA";
        } else if (s.equals("RC")) {
            return "Other";
        } else if (s.equals("repClass")) {
            return "Other";
        } else if (types.contains(s)) {
            return s;
        } else {
            return "Uncategorized";
        }

    }


}
