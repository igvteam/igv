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

package org.broad.igv.tools.converters;

import org.broad.igv.data.seg.Segment;

import java.io.*;
import java.util.*;

/**
 * User: jrobinso
 * Date: Apr 19, 2010
 */
public class IlluminaToCN {
    static int probeIdCol = 1;
    static int chrCol = 3;
    static int posCol = 4;
    static int preColCount = 10;
    static int sampleColumnCount = 16;
    static int cnvOffset = 5;
    static int logROffset = 7;
    static int afOffset = 8;

    /**
     * Converts an illumina aCGH file to multiple .seg and .cn files.  The files are produced
     *
     * @param iFile
     * @param oFilePrefix
     */

    public static void convertFile(String iFile, String oFilePrefix) {

        BufferedReader reader = null;
        PrintWriter cnvWriter = null;
        PrintWriter logRWriter = null;
        PrintWriter afWriter = null;

        LinkedHashMap<String, List<Segment>> cnvSegments = new LinkedHashMap();

        try {
            reader = new BufferedReader(new FileReader(iFile));
            cnvWriter = new PrintWriter(new BufferedWriter(new FileWriter(oFilePrefix + ".cnv.cn")));
            logRWriter = new PrintWriter(new BufferedWriter(new FileWriter(oFilePrefix + ".logr.cn")));
            afWriter = new PrintWriter(new BufferedWriter(new FileWriter(oFilePrefix + ".af.cn")));
            List<PrintWriter> allWriters = Arrays.asList(logRWriter, afWriter);

            String nextLine = reader.readLine();
            List<String> samples = getSamples(nextLine);

            for (PrintWriter pw : allWriters) {
                pw.print("Probe\tChr\tPosition");
                for (String s : samples) {
                    pw.print("\t" + s);
                }
                pw.println();
            }
            cnvWriter.println("Sample\tchr\tstart\tend\tnumberOfSnps");


            Map<String, Segment> currentSegments = new HashMap();

            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");

                String probe = tokens[probeIdCol];
                String chr = "chr" + tokens[chrCol];
                int position = Integer.parseInt(tokens[posCol]);


                for (PrintWriter pw : allWriters) {
                    pw.print(probe + "\t" + chr + "\t" + position);
                }

                int col = preColCount - 1;
                for (int i = 0; i < samples.size(); i++) {
                    String sample = samples.get(i);

                    int cnv = (int) Float.parseFloat(tokens[col + cnvOffset]);

                    Segment segment = currentSegments.get(sample);
                    if (segment == null || !chr.equals(segment.getChr()) || cnv == (int) segment.getScore()) {
                        if (segment != null) {
                            List<Segment> segs = cnvSegments.get(sample);
                            if (segs == null) {
                                segs = new ArrayList();
                                cnvSegments.put(sample, segs);
                            }
                            segs.add(segment);
                        }
                        segment = new Segment(position, position, cnv);
                        currentSegments.put(sample, segment);
                    } else {

                        //segment.incremenetSnpCount(1);
                        segment.setEnd(position);
                    }


                    String logR = tokens[col + logROffset];
                    logRWriter.print("\t" + logR);

                    String af = tokens[col + afOffset];
                    afWriter.print("\t" + af);

                    col += sampleColumnCount;
                }
                for (PrintWriter pw : allWriters) {
                    pw.println();
                }
            }

            for (Map.Entry<String, List<Segment>> entry : cnvSegments.entrySet()) {
                String sample = entry.getKey();
                for (Segment segment : entry.getValue()) {
                    cnvWriter.println(sample + "\t" + segment.getChr() + "\t" + segment.getStart() + "\t"
                            + segment.getEnd() + "\t" + "" + "\t" + (int) segment.getScore());
                }
            }

        }
        catch (IOException e) {
            e.printStackTrace();
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                }
            }
            for (PrintWriter pw : Arrays.asList(cnvWriter, logRWriter, afWriter)) {
                if (pw != null) {
                    pw.close();
                }
            }

        }

    }


    //SB516834_PB15681_E02.GType

    private static List<String> getSamples(String nextLine) {
        String[] tokens = nextLine.split("\t");
        int col = preColCount;

        List<String> samples = new ArrayList();
        while (col < tokens.length - 1) {
            samples.add(tokens[col].replace(".GType", ""));
            col += sampleColumnCount;
        }
        return samples;
    }

    public static void main(String[] args) {
        String ifile = "/Users/jrobinso/IGV/maggie/FullDataTable2.txt";
        String oFilePre = "/Users/jrobinso/IGV/maggie/smallDataTable";
        convertFile(ifile, oFilePre);
    }
}
