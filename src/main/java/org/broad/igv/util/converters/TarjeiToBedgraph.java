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
