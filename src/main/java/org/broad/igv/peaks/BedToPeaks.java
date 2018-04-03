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

package org.broad.igv.peaks;

import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.util.CompressionUtils;
import htsjdk.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.*;

/**
 * A utility class to convert a set of time-course bed files to a binary "peaks" file.
 *
 * @author jrobinso
 * @date Apr 22, 2011
 */
public class BedToPeaks {
    private static Map<String, String> colorMap;
    private static CompressionUtils compressionUtils;

    public BedToPeaks() {
        compressionUtils = new CompressionUtils();
    }

    /**
     * Converts a collection of peak "bed" files to a single binary "peak" file
     * <p/>
     * Note:  All input bed files must have identical structure.  Specifically
     * identical number of lines
     * for any line number, identical loci (chr, start, end)
     * identical placement of track lines and comments ("#" lines)
     */
    public static void createBinaryPeakFile(String factor, List<Integer> times, List<File> bedFiles, File outputDir, String root) throws IOException {

        String c = colorMap.get(factor);
        if (c == null) {
            System.out.println("No color found for " + factor);
            c = "0,0,150";
        }

        String peeksFileName = factor + ".peak.bin";
        File peeksFile = new File(outputDir, peeksFileName);


        BufferedReader[] readers = new BufferedReader[bedFiles.size()];
        LittleEndianOutputStream peakWriter = null;
        long indexPosition = 0l;
        try {

            for (int i = 0; i < bedFiles.size(); i++) {
                readers[i] = new BufferedReader(new FileReader(bedFiles.get(i)));
            }

            peakWriter = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(peeksFile)));

            LinkedHashMap<String, Long> chrIndex = new LinkedHashMap<String, Long>();

            String[] line = new String[readers.length];

            // Placeholder for index position
            peakWriter.writeLong(0l);

            // Track line
            peakWriter.writeString("track name=" + factor + " sample=" + factor + " viewLimits=0:100 useScore=1 color=" + c);

            // Time values
            int nTimePoints = times.size();
            peakWriter.writeInt(nTimePoints);
            for (int i = 0; i < times.size(); i++) {
                peakWriter.writeInt(times.get(i));
            }

            // Path to associated signal ("tdf") files
            String tdf = root + "/tdf/compressed/" + factor + ".merged.bam.tdf";
            peakWriter.writeString(tdf);
            for (int t : times) {
                peakWriter.writeString(root + "/tdf/timecourses/" + factor + "_" + t + "/" + factor + "_" + t + ".merged.bam.tdf");
            }


            // Process bed records
            String currentChr = "";
            List<PeakRecord> records = new ArrayList(20000);
            while (true) {

                for (int i = 0; i < readers.length; i++) {
                    line[i] = readers[i].readLine();
                    if (line[i] == null) {
                        // EOF, we're done
                        if (records.size() > 0) {
                            writeChromosomeData(peakWriter, chrIndex, currentChr, records);
                        }
                        break;
                    }
                }
                if (line[0].startsWith("#") || line[0].startsWith("track")) {
                    continue;   // Skip these lines.  Note we should actually check all files
                }

                // Take locus from first file
                String[] tokens = line[0].split("\t");
                String chr = tokens[0];


                if (!chr.equals(currentChr) && records.size() > 0) {
                    // End of data for current chromosome.
                    // Record position of start of this chromosome block for index
                    writeChromosomeData(peakWriter, chrIndex, currentChr, records);
                }

                if (chr == null) {
                    break;
                }

                currentChr = chr;

                int start = Integer.parseInt(tokens[1]);
                int end = Integer.parseInt(tokens[2]);
                float score = Float.parseFloat(tokens[4]);
                float[] timeScores = new float[nTimePoints];
                for (int i = 0; i < nTimePoints; i++) {
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

            // Write out index
            indexPosition = peakWriter.getWrittenCount();
            peakWriter.writeInt(chrIndex.size());
            for (Map.Entry<String, Long> entry : chrIndex.entrySet()) {
                peakWriter.writeString(entry.getKey());
                peakWriter.writeLong(entry.getValue());
            }

        } finally {
            if (peakWriter != null) {
                peakWriter.close();
                writeIndexPosition(peeksFile, indexPosition);
            }
            for (BufferedReader reader : readers) {
                reader.close();
            }
        }
    }

    private static void writeChromosomeData(LittleEndianOutputStream peakWriter,
                                            LinkedHashMap<String, Long> chrIndex,
                                            String currentChr,
                                            List<PeakRecord> records) throws IOException {

        chrIndex.put(currentChr, peakWriter.getWrittenCount());

        // Compress data for chromosome, then write it out.
        BufferedByteWriter buffer = new BufferedByteWriter(100000);

        buffer.putNullTerminatedString(currentChr);
        buffer.putInt(records.size());
        for (PeakRecord record : records) {
            buffer.putInt(record.start);
            buffer.putInt(record.end);
            buffer.putFloat(record.score);
            for (int i = 0; i < record.timeScores.length; i++) {
                buffer.putFloat(record.timeScores[i]);
            }
        }

        byte[] bytes = buffer.getBytes();
        bytes = compressionUtils.compress(bytes);
        peakWriter.writeInt(bytes.length);
        peakWriter.write(bytes);

        records.clear();
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

    static void writeIndexPosition(File file, long indexPosition) throws IOException {
        RandomAccessFile raf = null;
        try {
            raf = new RandomAccessFile(file, "rw");
            raf.getChannel().position(0);
            // Write as little endian
            BufferedByteWriter buffer = new BufferedByteWriter();
            buffer.putLong(indexPosition);
            raf.write(buffer.getBytes());
            raf.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            if (raf != null) raf.close();
        }
    }


    /**
     * Convenience method for processing files from the HT-Chip paper.
     *
     * @param rootDir
     * @param destDir
     * @param root
     * @throws IOException
     */
    public static void convertAll(File rootDir, File destDir, String root) throws IOException {
        // <FACTOR>.peaks.filtered.by.fold.real.sorted.bed
        // <FACTOR>_0.peaks.filtered.by.fold.real.sorted.bed

        int[] allTimes = {0, 30, 60, 120};

        File factorFile = new File(rootDir, "factors.txt");

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

            List<Integer> times = new ArrayList<Integer>();
            for (int i = 0; i < allTimes.length; i++) {
                bedFile = new File(rootDir, factor + "_" + allTimes[i] + ".peaks.filtered.by.fold.real.sorted.bed");
                if (bedFile.exists()) {
                    //System.out.println("Copying " + bedFile.getName());
                    bedFiles.add(bedFile);
                    times.add(allTimes[i]);
                } else {
                    System.out.println("Can't find " + bedFile.getName());
                }

            }

            createBinaryPeakFile(factor, times, bedFiles, destDir, root);
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


    /**
     * Process Broad ht-chip files.
     * <p/>
     * timePoints=0,30,60,120\
     *
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        File inputDir = new File(args[0] + "bed/");
        File destDir = new File(args[0] + args[1]);
        String root = args[2];
        colorMap = loadColors(args[0] + "colors.txt");
        convertAll(inputDir, destDir, root);
    }


}
