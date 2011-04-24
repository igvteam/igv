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

package org.broad.igv.peaks;

import org.broad.igv.tools.sort.Sorter;
import org.broad.igv.util.ColorUtilities;
import org.broad.igv.util.FileUtils;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A utitility class to convert a set of time-course bed files to a "peaks" file.  This is experimental.
 *
 * @author jrobinso
 * @date Apr 22, 2011
 */
public class BedToPeaks {
    private static Map<String, Color> colorMap;

    public static void main(String[] args) throws IOException {
        /*  List<File> bedFiles = new ArrayList();
       bedFiles.add(new File("/Users/jrobinso/IGV/time_course/cebpb.peaks.filtered.by.fold.real.sorted.bed"));
       bedFiles.add(new File("/Users/jrobinso/IGV/time_course/cebpb_0.peaks.filtered.by.fold.real.sorted.bed"));
       bedFiles.add(new File("/Users/jrobinso/IGV/time_course/cebpb_30.peaks.filtered.by.fold.real.sorted.bed"));
       bedFiles.add(new File("/Users/jrobinso/IGV/time_course/cebpb_60.peaks.filtered.by.fold.real.sorted.bed"));
       bedFiles.add(new File("/Users/jrobinso/IGV/time_course/cebpb_120.peaks.filtered.by.fold.real.sorted.bed"));

       File outputFile = new File("/Users/jrobinso/IGV/time_course/cebpb.peaks.igv");

       convert(bedFiles, outputFile);
        */

        //String rootDir = "/Volumes/seq_dcchip/mouse/DC/chipSeq/compressed/";
        //downloadAll(rootDir, new File("/Users/jrobinso/IGV/time_course"));
       // sortAll("/Users/jrobinso/IGV/time_course");

        File rootDir = new File("/Users/jrobinso/IGV/time_course");
        File inputDir = new File(rootDir, "sorted");
        File destDir = new File(rootDir, "peaks");

        convertAll(inputDir, destDir);
    }


    /**
     * Converts a collection of peak "bed" files to a single .igv file.
     *

     */
    public static void convert(String factor, List<File> bedFiles, File outputDir) throws IOException {


        if(colorMap == null) colorMap = loadColors("/Users/jrobinso/IGV/time_course/colors.txt");
        Color c = colorMap.get(factor);
        if(c == null) {
            System.out.println("No color found for " + factor);
            c = new Color(0, 0, 150);
        }


        BufferedReader[] readers = new BufferedReader[bedFiles.size()];
        PrintWriter writer = null;

        for (int i = 0; i < bedFiles.size(); i++) {
            readers[i] = new BufferedReader(new FileReader(bedFiles.get(i)));

        }

        File outputFile = new File(outputDir, factor + ".peak");
        writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

        writer.println("#track name=" + factor + " sample=" + factor + " viewLimits=0:100 useScore=1 color=" + ColorUtilities.colorToString(c));
        writer.print("Chr\tStart\tEnd\tName");
        for (File f : bedFiles) {
            writer.print("\t" + f.getName().replace(".bed", ""));
        }
        writer.println();


        // All input files must have the same # of lines, and be in the same order
        String[] line = new String[readers.length];

        try {
            while (true) {

                for (int i = 0; i < readers.length; i++) {
                    line[i] = readers[i].readLine();
                    if (line[i] == null) {
                        return;
                    }
                }

                if (line[0].startsWith("#") || line[0].startsWith("track")) {
                    continue;
                }

                // Take locus from first file
                String[] tokens = line[0].split("\t");
                String chr = tokens[0];
                int start = Integer.parseInt(tokens[1]);
                int end = Integer.parseInt(tokens[2]);
                String name = tokens[3];
                float score = Float.parseFloat(tokens[4]);
                writer.print(chr + "\t" + start + "\t" + end + "\t" + name + "\t" + score);
                for (int i = 1; i < line.length; i++) {
                    tokens = line[i].split("\t");
                    if (!(tokens[0].equals(chr) &&
                            Integer.parseInt(tokens[1]) == start &&
                            Integer.parseInt(tokens[2]) == end)) {
                        throw new RuntimeException("Unordered files");
                    }
                    score = Float.parseFloat(tokens[4]);
                    writer.print("\t" + score);
                }
                writer.println();
            }
        } finally {
            writer.close();
            for (BufferedReader reader : readers) {
                reader.close();
            }
        }
    }


    public static void sortAll(String rootDir) {

        for (File f : (new File(rootDir)).listFiles()) {
            if (f.getName().endsWith(".bed") && !f.getName().endsWith(".sorted.bed")) {
                File of = new File(f.getAbsolutePath().replace(".bed", ".sorted.bed"));
                if (!of.exists()) {
                    System.out.println("Sorting " + f.getName());
                    Sorter sorter = Sorter.getSorter(f, of);
                    sorter.run();
                }
            }
        }

    }

    public static void convertAll(File rootDir, File destDir) throws IOException {
        // <FACTOR>.peaks.filtered.by.fold.real.sorted.bed
        // <FACTOR>_0.peaks.filtered.by.fold.real.sorted.bed


        File factorFile = new File("/Users/jrobinso/IGV/time_course/factors.txt");
        int[] time = {0, 30, 60, 120};

        BufferedReader reader = new BufferedReader(new FileReader(factorFile));
        String factor;
        while ((factor = reader.readLine()) != null) {
            factor = factor.trim();

            List<File> bedFiles = new ArrayList();

            File bedFile = new File(rootDir, factor + ".peaks.filtered.by.fold.real.sorted.bed");
            if (bedFile.exists()) {
                //System.out.println("Adding " + bedFile.getName());
                bedFiles.add(bedFile);
            }
            else {
                System.out.println("Can't find " + bedFile.getName());
            }


            for (int i = 0; i < time.length; i++) {
                bedFile = new File(rootDir, factor + "_" + time[i] + ".peaks.filtered.by.fold.real.sorted.bed");
                if (bedFile.exists()) {
                    //System.out.println("Copying " + bedFile.getName());
                    bedFiles.add(bedFile);
                }
                else {
                    System.out.println("Can't find " + bedFile.getName());
                }

            }

            convert(factor, bedFiles, destDir);


        }
        reader.close();


    }

    public static void downloadAll(String rootDir, File destDir) throws IOException {
        // <FACTOR>/segments_0.05/<FACTOR>.peaks.filtered.by.fold.real.bed
        // timecourses/<FACTOR>_0/segments_0.05/<FACTOR>_0.peaks.filtered.by.fold.real.bed

        File factorFile = new File("/Users/jrobinso/IGV/time_course/factors.txt");
        int[] time = {0, 30, 60, 120};

        BufferedReader reader = new BufferedReader(new FileReader(factorFile));
        String factor;
        while ((factor = reader.readLine()) != null) {
            factor = factor.trim();

            List<File> bedFiles = new ArrayList();

            File bedFile = new File(rootDir + factor + "/segments_0.05/" + factor + ".peaks.filtered.by.fold.real.bed");
            if (bedFile.exists()) {
                System.out.println("Copying " + bedFile.getName());
                FileUtils.copyFile(bedFile, new File(destDir, bedFile.getName()));
            }


            for (int i = 0; i < time.length; i++) {
                bedFile = new File(rootDir + "timecourses/" + factor + "_" + time[i] + "/segments_0.05/" +
                        factor + "_" + time[i] + ".peaks.filtered.by.fold.real.bed");
                if (bedFile.exists()) {
                    System.out.println("Copying " + bedFile.getName());
                    FileUtils.copyFile(bedFile, new File(destDir, bedFile.getName()));
                }

            }

        }
        reader.close();
    }

    public static Map<String, Color> loadColors(String file) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        Map<String, Color> map = new HashMap();
        while ((line = reader.readLine()) != null) {
            String [] tokens = line.split("\t");
            Color c = ColorUtilities.stringToColor(tokens[1]);
            map.put(tokens[0], c);
        }
        reader.close();
        return map;
    }

}
