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

package org.broad.igv.peaks;

import org.broad.igv.tools.sort.Sorter;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

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
    private static Map<String, String> colorMap;

    public static void main(String[] args) throws IOException {

        File inputDir = new File("/Users/jrobinso/IGV/ichip/bed/");
        File destDir = new File("/Users/jrobinso/IGV/ichip/peaks");

       // downloadAll("/Volumes/seq_dcchip/mouse/DC/chipSeq/compressed/", inputDir);
       // sortAll("/Users/jrobinso/IGV/time_course/sorted");
        convertAll(inputDir, destDir);
    }


    /**
     * Converts a collection of peak "bed" files to a single .igv file.
     */
    public static void createCfgFile(String factor, File outputDir) throws IOException {

        int[] times = {0, 30, 60, 120};

        if (colorMap == null) colorMap = loadColors("/Users/jrobinso/IGV/time_course/colors.txt");
        String c = colorMap.get(factor);
        if (c == null) {
            System.out.println("No color found for " + factor);
            c = "0,0,150";
        }

        PrintWriter cfgWriter = null;
        try {

            File cfgFile = new File(outputDir, factor + ".peak.cfg");
            cfgWriter = new PrintWriter(new BufferedWriter(new FileWriter(cfgFile)));

            String peeksFile = factor + ".peak";

            // Signals at http://iwww.broadinstitute.org/igvdata/peaks/compressed/<FACTOR>/<FACTOR>.merged.bam.tdf
            String peaks = "http://www.broadinstitute.org/igvdata/ichip/peaks/" + peeksFile;
            String tdf = "http://www.broadinstitute.org/igvdata/ichip/tdf/compressed/" + factor + ".merged.bam.tdf";
            String bam = "http://www.broadinstitute.org/igvdata/ichip/compressed/" + factor + "/" + factor + ".merged.bam";

            cfgWriter.println("track name=" + factor + " sample=" + factor + " viewLimits=0:100 useScore=1 color=" + c);
            cfgWriter.println("timePoints=0,30,60,120");
            cfgWriter.println("peaks=" + peaks);
            cfgWriter.println("signals=" + tdf);
            cfgWriter.print("timeSignals=");

            String root = "http://www.broadinstitute.org/igvdata/ichip/tdf/timecourses/";
            for (int t : times) {
                cfgWriter.print(root + factor + "_" + t + "/" + factor + "_" + t + ".merged.bam.tdf,");
            }
            cfgWriter.println();
            cfgWriter.println("#index");

            File pf = new File(outputDir, peeksFile);
            indexPeakFile(pf, cfgWriter);


        } finally {
            if (cfgWriter != null) cfgWriter.close();
        }
    }

    private static void indexPeakFile(File peeksFile, PrintWriter cfgWriter) {
        AsciiLineReader reader = null;

        try {
            reader = ParsingUtils.openAsciiReader(new ResourceLocator(peeksFile.getAbsolutePath()));
            String nextLine;
            // Skip first 2 lines (track line and header line)
            reader.readLine();
            reader.readLine();
            String lastChr = "";
            long position = reader.getPosition();
            while ((nextLine = reader.readLine()) != null) {
                String chr = nextLine.split("\t")[0];
                if (!chr.equals(lastChr)) {
                    cfgWriter.println(chr + "\t" + position);
                }
                position = reader.getPosition();
                lastChr = chr;
            }


        }
        catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } finally {

        }
    }


    /**
     * Converts a collection of peak "bed" files to a single .igv file.
     */
    public static void createBinaryPeakFile(String factor, List<File> bedFiles, File outputDir) throws IOException {


        BufferedReader[] readers = new BufferedReader[bedFiles.size()];
        DataOutputStream peakWriter = null;
        try {

            for (int i = 0; i < bedFiles.size(); i++) {
                readers[i] = new BufferedReader(new FileReader(bedFiles.get(i)));
            }

            File peeksFile = new File(outputDir, factor + ".peak.bin");
            peakWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(peeksFile)));

            if (peakWriter != null) {


                // All input files must have the same # of lines, and be in the same order
                String[] line = new String[readers.length];

                int nTimePoints = readers.length - 1;
                peakWriter.writeInt(nTimePoints);
                String lastChr = "";
                List<PeakRecord> records = new ArrayList(20000);
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

                    if (!chr.equals(lastChr) && records.size() > 0) {
                        peakWriter.writeUTF(lastChr);
                        peakWriter.writeInt(records.size());
                        for (PeakRecord record : records) {
                            peakWriter.writeInt(record.start);
                            peakWriter.writeInt(record.end);
                            peakWriter.writeFloat(record.score);
                            for (int i = 0; i < record.timeScores.length; i++) {
                                peakWriter.writeFloat(record.timeScores[i]);
                            }
                        }
                        records.clear();
                    }
                    lastChr = chr;

                    int start = Integer.parseInt(tokens[1]);
                    int end = Integer.parseInt(tokens[2]);
                    float score = Float.parseFloat(tokens[4]);
                    float [] timeScores = new float[nTimePoints];
                    for (int i = 0; i<nTimePoints; i++) {
                        tokens = line[i + 1].split("\t");
                        if (!(tokens[0].equals(chr) &&
                                Integer.parseInt(tokens[1]) == start &&
                                Integer.parseInt(tokens[2]) == end)) {
                            throw new RuntimeException("Unordered files");
                        }
                        timeScores[i] = Float.parseFloat(tokens[4]);
                    }
                    records.add(new PeakRecord(start, end, score, timeScores));
                }
            }

        } finally {
            if (peakWriter != null) {
                peakWriter.writeUTF("EOF");
                peakWriter.close();
            }
            for (BufferedReader reader : readers) {
                reader.close();
            }
        }
    }

    static class PeakRecord {
        int start;
        int end;
        float score;
        float[] timeScores;

        PeakRecord(int start, int end, float score, float[] timeScores) {
            this.start = start;
            this.end = end;
            this.score = score;
            this.timeScores = timeScores;
        }
    }


    /**
     * Converts a collection of peak "bed" files to a single .igv file.
     */
    public static void createPeakFile(String factor, List<File> bedFiles, File outputDir) throws IOException {

        int[] times = {0, 30, 60, 120};

        if (colorMap == null) colorMap = loadColors("/Users/jrobinso/IGV/time_course/colors.txt");
        String c = colorMap.get(factor);
        if (c == null) {
            System.out.println("No color found for " + factor);
            c = "0,0,150";
        }

        BufferedReader[] readers = new BufferedReader[bedFiles.size()];
        PrintWriter peakWriter = null;
        PrintWriter cfgWriter = null;
        try {

            for (int i = 0; i < bedFiles.size(); i++) {
                readers[i] = new BufferedReader(new FileReader(bedFiles.get(i)));

            }

            File peeksFile = new File(outputDir, factor + ".peak");
            peakWriter = new PrintWriter(new BufferedWriter(new FileWriter(peeksFile)));


            // Signals at http://iwww.broadinstitute.org/igvdata/peaks/compressed/<FACTOR>/<FACTOR>.merged.bam.tdf
            String peaks = "http://www.broadinstitute.org/igvdata/ichip/peaks/" + peeksFile.getName();
            String tdf = "http://www.broadinstitute.org/igvdata/ichip/tdf/compressed/" + factor + ".merged.bam.tdf";
            String bam = "http://www.broadinstitute.org/igvdata/ichip/compressed/" + factor + "/" + factor + ".merged.bam";
            //String localBam = "/seq/dcchip/mouse/DC/chipSeq/compressed/"+ factor + "/" + factor + ".merged.bam";
            //System.out.println("bsub -o " + factor + "out.txt  -q genepattern  \"/xchip/igv/tools/Test/igvtools count -w 1 -x t=2 " +
            //       "-c " + c + " " + localBam + " " + factor + ".tdf mm9\"");
            //System.out.println("        <Resource name=\"" + factor + "\"  trackLine=' + trackLine + "'\n" +
            //         "                  path=\"" + tdf + "\"/>");

            //if(true) return;


            if (peakWriter != null) {
                peakWriter.println("track name=" + factor + " sample=" + factor + " viewLimits=0:100 useScore=1 color=" + c);

                //peakWriter.print("Chr\tStart\tEnd\tName");
                for (File f : bedFiles) {
                    peakWriter.print("\t" + f.getName().replace(".bed", ""));
                }
                peakWriter.println();


                // All input files must have the same # of lines, and be in the same order
                String[] line = new String[readers.length];

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
                    peakWriter.print(chr + "\t" + start + "\t" + end + "\t" + name + "\t" + score);
                    for (int i = 1; i < line.length; i++) {
                        tokens = line[i].split("\t");
                        if (!(tokens[0].equals(chr) &&
                                Integer.parseInt(tokens[1]) == start &&
                                Integer.parseInt(tokens[2]) == end)) {
                            throw new RuntimeException("Unordered files");
                        }
                        score = Float.parseFloat(tokens[4]);
                        peakWriter.print("\t" + score);
                    }
                    peakWriter.println();
                }
            }

        } finally {
            if (peakWriter != null) peakWriter.close();
            if (cfgWriter != null) cfgWriter.close();
            for (BufferedReader reader : readers) {
                reader.close();
            }
        }
    }


    public static void sortAll(String rootDir) throws IOException {

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


        File factorFile = new File("/Users/jrobinso/IGV/ichip/bed/factors.txt");
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
            } else {
                System.out.println("Can't find " + bedFile.getName());
            }


            for (int i = 0; i < time.length; i++) {
                bedFile = new File(rootDir, factor + "_" + time[i] + ".peaks.filtered.by.fold.real.sorted.bed");
                if (bedFile.exists()) {
                    //System.out.println("Copying " + bedFile.getName());
                    bedFiles.add(bedFile);
                } else {
                    System.out.println("Can't find " + bedFile.getName());
                }

            }

            createBinaryPeakFile(factor, bedFiles, destDir);
            //createCfgFile(factor, destDir);


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

            //_0.01 rather than the segments_0.05 for chromatin (K4me3, K4me1 & K27Ac).
            String segments;
            if(factor.equals("K4me3") || factor.equals("K4me1") || factor.equals("K27Ac")) {
                 segments = "/segments_0.01/";
            }
            else {
                segments = "/segments_0.05/";
            }

            File bedFile = new File(rootDir + factor + segments + factor + ".peaks.filtered.by.fold.real.bed");
            if (bedFile.exists()) {
                System.out.println("Copying " + bedFile.getName());
                FileUtils.copyFile(bedFile, new File(destDir, bedFile.getName()));
            }
            else {
                System.out.println("File not found: " + bedFile);
            }


            for (int i = 0; i < time.length; i++) {
                bedFile = new File(rootDir + "timecourses/" + factor + "_" + time[i] + segments +
                        factor + "_" + time[i] + ".peaks.filtered.by.fold.real.bed");
                if (bedFile.exists()) {
                    System.out.println("Copying " + bedFile.getName());
                    FileUtils.copyFile(bedFile, new File(destDir, bedFile.getName()));
                }

            }

        }
        reader.close();
    }

    public static Map<String, String> loadColors(String file) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        Map<String, String> map = new HashMap();
        while ((line = reader.readLine()) != null) {
            String[] tokens = line.split("\t");
            map.put(tokens[0], tokens[1]);
        }
        reader.close();
        return map;
    }

}
