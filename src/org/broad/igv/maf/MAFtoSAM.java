/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
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

package org.broad.igv.maf;

import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
import java.io.*;


/**
 * Created by jrobinson on 5/5/16.
 * <p/>
 * Parses a MAF file for Alignment records, as opposed to Multiple Alignments.  For LAST aligner output.
 */

public class MAFtoSAM {


    public static void main(String[] args) throws IOException {
        String inputPath = args[0];
        String outputPath = args.length > 1 ? args[1] : null;
        convert(inputPath, outputPath);
    }


    public static void convert(String path, String outputPath) throws IOException {


        BufferedReader reader = ParsingUtils.openBufferedReader(path);

        PrintWriter out = outputPath == null ? new PrintWriter(System.out) :
                new PrintWriter(new BufferedWriter(new FileWriter(outputPath)));

        String line;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("a ")) {
                // Parse alignments until blank line
                parseBlock(reader, out);
            }
        }

        out.flush();
        out.close();
    }

    private static void parseBlock(BufferedReader reader, PrintWriter out) throws IOException {

        String line;
        SequenceLine referenceLine = null;
        SequenceLine queryLine;
        byte[] refBytes = null;

        while ((line = reader.readLine()) != null) {
            if (line.trim().length() == 0) {
                return; // return someething
            }
            if (line.startsWith("s ")) {
                if (null == referenceLine) {
                    referenceLine = parseSequenceLine(line);
                    refBytes = referenceLine.text.getBytes();
                } else {
                    queryLine = parseSequenceLine(line);

                    // Build cigar string one byte at a time
                    byte[] queryBytes = queryLine.text.getBytes();

                    if (queryBytes.length != refBytes.length)
                        throw new RuntimeException("Query and ref bytes unequal length");

                    String cigarString = "";

                    for (int i = 0; i < queryBytes.length; i++) {
                        byte q = queryBytes[i];
                        byte ref = refBytes[i];

                        if (q == '-') {
                            if (ref == '-') {
                                //ignore, caused by insertion in another alignment;
                            } else {
                                cigarString += "D";
                            }
                        } else {
                            if (ref == '-') {
                                cigarString += "I";
                            } else {
                                cigarString += "M";
                            }
                        }
                    }

                    cigarString = collapseCigar(cigarString);

                    // Output SAM record

                    String chr;
                    if(referenceLine.src.contains(".")) {
                        int idx = referenceLine.src.lastIndexOf('.') + 1;
                        chr = referenceLine.src.substring(idx);
                    }
                    else {
                        chr = referenceLine.src;
                    }

                    out.println(queryLine.src + "\t" +
                            0 + "\t" +
                            chr + "\t" +
                            (referenceLine.start + 1) + "\t" +
                            30 + "\t" +
                            cigarString + "\t" +
                            "*\t0\t0\t*\t*");
                }
            }
        }
    }

    private static String collapseCigar(String cigarString) {

        if (cigarString.length() == 0) return "";

        String collapsedCigar = "";
        char lastOperator = cigarString.charAt(0);
        int counter = 1;
        for (int i = 1; i < cigarString.length(); i++) {
            if (cigarString.charAt(i) == lastOperator) {
                counter++;
            } else {
                collapsedCigar += ("" + counter + lastOperator);
                lastOperator = cigarString.charAt(i);
                counter = 1;
            }
        }
        collapsedCigar += ("" + counter + lastOperator);
        return collapsedCigar;
    }


    /**
     * Parse an alignment block.   Assumes 1 alignment per block, first "s" record is reference, second is alignment
     *
     * @param line
     */

    private static SequenceLine parseSequenceLine(String line) throws IOException {


        String[] tokens = Globals.whitespacePattern.split(line);
        SequenceLine sl = new SequenceLine();
        sl.src = tokens[1];
        sl.start = Integer.parseInt(tokens[2]);
        sl.size = Integer.parseInt(tokens[3]);
        sl.strand = tokens[4].charAt(0);
        sl.srcSize = Integer.parseInt(tokens[5]);
        sl.text = tokens[6];
        return sl;

    }

    static class SequenceLine {
        String src;
        int start;
        int size;
        char strand;
        int srcSize;
        String text;
    }
}

