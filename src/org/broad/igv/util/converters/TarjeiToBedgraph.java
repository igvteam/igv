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

package org.broad.igv.util.converters;

import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;

/**
 * Converts a director of files of the form chr1.Brain.H3K4me2.txt in "Tarjei" format to bedgraph format.
 *
 * @author Jim Robinson
 * @date 9/13/11
 */
public class TarjeiToBedgraph {

    public static void main(String[] args) throws IOException {

        File dir = new File(args[0]);
        File[] files = dir.listFiles();

        // Parse out chromosome and sample names
        LinkedHashSet<String> chromosomes = new LinkedHashSet<String>();
        LinkedHashSet<String> samples = new LinkedHashSet<String>();

        for (File f : files) {
            String[] tokens = f.getName().split("\\.");
            chromosomes.add(tokens[0]);
            String sample = tokens[1];
            for (int idx = 2; idx < tokens.length - 1; idx++) {
                sample += ("." + tokens[idx]);
            }
            samples.add(sample);
        }

        for (String sample : samples) {
            LinkedHashMap<String, String> chrFileMap = new LinkedHashMap<String, String>();
            for (String chr : chromosomes) {
                File f = new File(dir, chr + "." + sample + ".txt");
                if (f.exists()) {
                    chrFileMap.put(chr, f.getAbsolutePath());
                }
            }

            String ofile = (new File(dir, sample + ".bedgraph")).getAbsolutePath();
            convert(chrFileMap, ofile);
        }
    }

    public static void convert(LinkedHashMap<String, String> chrFileMap, String ofile) throws IOException {
        int step = 25;

        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(ofile)));

        for (Map.Entry<String, String> entry : chrFileMap.entrySet()) {

            String chr = entry.getKey();
            String ifile = entry.getValue();
            BufferedReader reader = null;

            reader = ParsingUtils.openBufferedReader(ifile);

            String nextLine;
            int start = 0;
            while ((nextLine = reader.readLine()) != null) {
                float counts = Float.parseFloat(nextLine.trim());
                if (counts > 0) {
                    pw.println(chr + "\t" + start + "\t" + (start + step) + "\t" + counts);
                }
                start += step;
            }

            reader.close();
        }
        pw.close();
    }

}
