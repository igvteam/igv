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

package org.broad.igv.feature;


import org.broad.igv.Globals;


import java.io.*;
import java.lang.reflect.Array;
import java.util.*;


class LocAndVal {
    public int loc;
    public double val;
    public LocAndVal(int loc, double val) {
        this.loc = loc;
        this.val = val;
    }
}



/**
 * @author sbusan
 */
public class ShapeFileUtils {

    static LinkedList<LocAndVal> transformProfile(LinkedList<LocAndVal> profile,
                                                  int seqLen,
                                                  int newLeft,
                                                  String strand){
        LinkedList<LocAndVal> transProfile = new LinkedList<>();
        for (LocAndVal d : profile){
            int loc = d.loc;
            double val = d.val;
            if (strand == "+"){
                loc = loc + newLeft - 1;
            } else if (strand == "-"){
                loc = seqLen - loc + newLeft;
            } else {
                throw new RuntimeException("Unrecognized strand (options: \"+\",\"-\")");
            }
            transProfile.add(new LocAndVal(loc, val));
        }
        return transProfile;
    }

    static LinkedList<LocAndVal> loadShape(String inFile) throws
            FileNotFoundException, IOException {
        // TODO: add error messages for misformatted file
        LinkedList<LocAndVal> profile = new LinkedList<>();

        BufferedReader br = null;

        try {
            br = new BufferedReader(new FileReader(inFile));

            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                String[] s = Globals.whitespacePattern.split(nextLine.trim());
                int loc = Integer.parseInt(s[0]);
                double val = Double.parseDouble(s[1]);
                if (val < -998) {
                    val = Double.NaN;
                }
                profile.add(new LocAndVal(loc, val));
            }
        } finally {
            if (br != null) br.close();
        }

        return profile;
    }


    static void writeWigFile(String wigFile,
                             String chrom,
                             LinkedList<LocAndVal> profile) throws IOException {
        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(wigFile)));
            // write header
            pw.println("variableStep chrom="+chrom);

            // write locs and values
            for (LocAndVal d : profile) {
                int loc = d.loc;
                double val = d.val;
                if (!Double.isNaN(val)){
                    pw.println(""+loc+"\t"+String.format("%.6f",val));
                }
            }
        } finally {
            if (pw != null) pw.close();
        }
    }


    public static void shapeToWigFile(String inFile,
                                      String wigFile,
                                      String chromosome,
                                      String strand,
                                      int left) throws
            FileNotFoundException, IOException {

        LinkedList<LocAndVal> profile = loadShape(inFile);
        profile = transformProfile(profile, profile.size(), left, strand);
        writeWigFile(wigFile, chromosome, profile);
    }

    /**
     * Convert a base pairing structure file in dot-bracket notation
     * (also known as Vienna format) to an easily parseable .bp arcs file. Does not
     * currently handle mapping coords to spliced transcripts.
     *
     * @param dbFile        Input file
     * @param bpFile        Output file
     * @param chromosome    Associated IGV chromosome
     * @param strand        Associated strand ("+" or "-")
     * @param left          Starting left-most position (0-based)
     */


    /**
     * Convert a pairing probability file as output by RNAStructure
     * and/or SuperFold to an easily parseable .bp arcs file. Does not
     * currently handle mapping coords to spliced transcripts.
     *
     * @param dpFile        Input file
     * @param bpFile        Output file
     * @param chromosome    Associated IGV chromosome
     * @param strand        Associated strand ("+" or "-")
     * @param left          Starting left-most position (0-based)
     */

}
