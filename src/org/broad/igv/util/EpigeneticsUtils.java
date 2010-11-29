/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import org.broad.igv.track.TrackType;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * Utility methods to deal with Bernstein's groups file conventions.
 */

public class EpigeneticsUtils {


    static int step = 25;
    static int span = 25;

    public static void main(String[] args) {
        try {
            String dir = "/Volumes/igv/data/public/epigenetics/WilmsTumor/wig/orig";
            String outputDir = "/Volumes/igv/data/public/epigenetics/WilmsTumor/wig";
            //convertToTDF(new File(dir), outputDir);
            createCombinedWigs(new File(dir), new File(outputDir));
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    public static void convertToTDF(File inputDir, String outputDir) {
        for (File f : inputDir.listFiles()) {
            String ifile = f.getName();
            if (ifile.endsWith(".wig")) {
                String oFile = (outputDir + f.getName() + ".tdf");
                String logFile = f.getName() + ".out";
                System.out.println("/xchip/igv/tools/igvtools tile " + ifile + " " + oFile + " hg18 > " + logFile + " &");
            }
        }
    }



    public static void createCombinedWigs(File inputDir, File outputDir) throws IOException {


        Map<String, List<File>> fileMap = getFilesBySample(inputDir);

        for (Map.Entry<String, List<File>> entry : fileMap.entrySet()) {

            BufferedReader br;
            PrintWriter pw;

            String sampleName = entry.getKey();
            File oFile = new File(outputDir, sampleName + ".wig");

            pw = new PrintWriter(new BufferedWriter(new FileWriter(oFile)));

            String cs = getColor(sampleName);
            pw.println("#type=" + TrackType.CHIP);
            pw.print("track  type=wiggle_0 name=" + sampleName + " viewLimits=0:10");
            if (cs == null) {
                pw.println();
            } else {
                pw.println(" color=" + cs);
            }

            for (File f : entry.getValue()) {

                br = ParsingUtils.openBufferedReader(f.getAbsolutePath());
                String stepLine = br.readLine();
                if (!stepLine.startsWith("fixedStep") || !stepLine.contains("start=1")) {
                    throw new RuntimeException("Unexpected step line: " + stepLine);
                }

                String[] tokens = stepLine.split(" ");
                String chromString = null;
                for (String s : tokens) {
                    if (s.startsWith("chrom")) {
                        chromString = s.trim();
                        break;
                    }
                }
                if (chromString == null) {
                    throw new RuntimeException("Missing chrom string: " + stepLine);
                }
                String chr = chromString.split("=")[1];





                // First line
                String nextLine = br.readLine();
                int value = Integer.parseInt(nextLine.trim());
                int startPos = 0;
                int endPos = step;

                while ((nextLine = br.readLine()) != null) {

                    int v = Integer.parseInt(nextLine.trim());

                    if (v != value && value >= 0) {
                        pw.println(chr + "\t" + startPos + "\t" + endPos + "\t" + value);
                        startPos = endPos;

                    }

                    if (value < 0) {
                        startPos += step;
                    }
                    value = v;
                    endPos += step;
                }

                // Write the last interval
                if (value >= 0) {
                    pw.println(chr + "\t" + startPos + "\t" + endPos + "\t" + value);
                    startPos = endPos;

                }


                br.close();

                /*
                pw.println("variableStep " + chromString + " span=" + span);
                int pos = 1;
                String nextLine;
                while ((nextLine = br.readLine()) != null) {
                    int value = Integer.parseInt(nextLine.trim());
                    if (value > 0) {
                        pw.println(pos + "\t" + value);
                    }
                    pos += step;
                }
                br.close();
                */
            }


            pw.close();


        }

    }


    /**
     * @param directory
     * @return
     */
    static Map<String, List<File>> getFilesBySample(File directory) {

        Map<String, List<File>> fileMap = new HashMap();

        for (File f : directory.listFiles()) {

            String fn = f.getName();
            if (fn.endsWith("wig") || fn.endsWith("wig.gz")) {
                String[] tokens = fn.split("_");
                if (tokens.length > 0) {
                    int sz = tokens.length - 1;
                    String s = tokens[0];
                    for (int i = 1; i < sz; i++) {
                        s += "_" + tokens[i];
                    }

                    List<File> files = fileMap.get(s);
                    if (files == null) {
                        files = new ArrayList();
                        fileMap.put(s, files);
                    }
                    files.add(f);

                }
            }

        }

        return fileMap;
    }


    static String getColor(String sampleName) {

        if (sampleName.contains("K4")) {
            return "0,150,0";
        } else if (sampleName.contains("K9")) {
            return "100,0,0";
        } else if (sampleName.contains("K27")) {
            return "255,0,0";
        }
        return null;


    }


}
