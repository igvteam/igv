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

import htsjdk.samtools.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.ChromosomeNameComparator;
import org.broad.igv.tools.sort.SAMSorter;
import org.broad.igv.tools.sort.Sorter;
import org.broad.igv.tools.sort.SorterFactory;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.*;


/**
 * Created by jrobinson on 5/5/16.
 * <p>
 * Parses a MAF file for Alignment records, as opposed to Multiple Alignments.  For LAST aligner output.
 */

public class MAFtoSAM {

    public static final int READ_PAIRED_FLAG = 0x1;
    public static final int PROPER_PAIR_FLAG = 0x2;
    public static final int READ_UNMAPPED_FLAG = 0x4;
    public static final int MATE_UNMAPPED_FLAG = 0x8;
    public static final int READ_STRAND_FLAG = 0x10;
    public static final int MATE_STRAND_FLAG = 0x20;
    public static final int FIRST_OF_PAIR_FLAG = 0x40;
    public static final int SECOND_OF_PAIR_FLAG = 0x80;
    public static final int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    public static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    public static final int DUPLICATE_READ_FLAG = 0x400;
    public static final int SUPPLEMENTARY_FLAG = 0x800;


    public static void main(String[] args) throws IOException {
        String inputPath = args[0];
        String outputPath = args.length > 1 ? args[1] : args[0] + ".sam";
        convert(inputPath, outputPath, false);
    }


    public static void convert(String path, String outputPath, boolean noSATag) throws IOException {

        File unsortedOutput = new File(outputPath + ".unsorted.sam");
        File outputFile = new File(outputPath);

        BufferedReader reader = ParsingUtils.openBufferedReader(path);
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(unsortedOutput)));
        Map<String, Integer> sequenceDictionary = new LinkedHashMap<String, Integer>();
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("a ")) {
                // Parse alignments until blank line
                parseBlock(reader, out, sequenceDictionary);
            }
        }
        out.flush();
        out.close();
        reader.close();

        // Insert header and  SA tags
        addSATags(unsortedOutput, outputFile, sequenceDictionary);


       unsortedOutput.deleteOnExit();

    }

    private static void addSATags(File inputFile, File outputFile, Map<String, Integer> sequenceDictionary) throws IOException {

        String groupedOutput = inputFile.getAbsolutePath() + ".grouped.sam";
        Sorter sorter = SorterFactory.getSorter(inputFile, new File(groupedOutput));
        sorter.setComparator(SAMSorter.ReadNameComparator);
        sorter.run();

        // Get Sam file reader, output header, then iterate through sam records combining those with like read names

        BufferedReader reader = new BufferedReader(new FileReader(groupedOutput));
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
        outputHeader(sequenceDictionary, out);


        String line;
        String[] lastRecord = null;

        List<String[]> chimericAlignments = new ArrayList<>();

        while ((line = reader.readLine()) != null) {
            if (line.startsWith("@")) {
                out.println(line);
            } else {
                String[] record = Globals.tabPattern.split(line);

                if (lastRecord == null || lastRecord[0].equals(record[0])) {
                    chimericAlignments.add(record);
                } else {
                    if (chimericAlignments.size() > 0) {
                        printChimericAlignments(out, chimericAlignments);
                    }
                    chimericAlignments.clear();
                    chimericAlignments.add(record);
                }
                lastRecord = record;
            }
        }

        // Last one
        if (chimericAlignments.size() > 0) {
            printChimericAlignments(out, chimericAlignments);
        }

        out.flush();
        out.close();

        (new File(groupedOutput)).deleteOnExit();

    }

    private static void printChimericAlignments(PrintWriter out, List<String[]> chimericAlignments) {

        if (chimericAlignments.size() == 1) {
            // Not chimeric, just print
            printRecord(out, chimericAlignments.get(0));
            out.println();
        } else {

            for (int i = 0; i < chimericAlignments.size(); i++) {

                String[] thisRecord = chimericAlignments.get(i);

                StringBuffer saTag = new StringBuffer("SA:Z:");


                // Build up the SA string -- essentially a list of the other parts of this chmeric alignment.
                for (int j = 0; j < chimericAlignments.size(); j++) {

                    if (j == i) continue;

                    String[] supRecord = chimericAlignments.get(j);

                    // This is fragile, but we control the order
                    String nm;
                    String nmString = supRecord[11];
                    if (nmString.startsWith("NM:i:")) {
                        nm = nmString.substring(5);
                    } else {
                        nm = "0";
                    }

                    boolean isNegativeStrand = (Integer.parseInt(supRecord[1]) & READ_STRAND_FLAG) != 0;

                    saTag.append(supRecord[2] + "," + supRecord[3] + "," + (isNegativeStrand ? "-," : "+,") +
                            supRecord[5] + "," + supRecord[4] + "," + nm + ";");
                }

                if (i > 0) {
                    int flags = Integer.parseInt(thisRecord[1]) | SUPPLEMENTARY_FLAG;
                    thisRecord[1] = Integer.toString(flags);
                }

                printRecord(out, thisRecord);

                out.print("\t" +saTag.toString());

                out.println();


            }
        }

    }

    private static void printRecord(PrintWriter out, String[] record) {
        out.print(record[0]);
        for (int i = 1; i < record.length; i++) {
            out.print("\t" + record[i]);
        }
    }


    private static void parseBlock(BufferedReader reader, PrintWriter out, Map<String, Integer> sequenceDictionary) throws IOException {

        String line;
        SequenceLine referenceLine = null;
        SequenceLine queryLine;
        byte[] refBytes = null;
        String chr = null;

        while ((line = reader.readLine()) != null) {
            if (line.trim().length() == 0) {
                return; // return someething
            }
            if (line.startsWith("s ")) {

                if (null == referenceLine) {
                    referenceLine = parseSequenceLine(line);
                    refBytes = referenceLine.text.getBytes();

                    if (referenceLine.src.contains(".")) {
                        int idx = referenceLine.src.lastIndexOf('.') + 1;
                        chr = referenceLine.src.substring(idx);
                    } else {
                        chr = referenceLine.src;
                    }
                    sequenceDictionary.put(chr, referenceLine.srcSize);
                } else {

                    queryLine = parseSequenceLine(line);

                    int flags = queryLine.strand == '+' ? 0 : 16;

                    // Build cigar string one byte at a time
                    byte[] queryBytes = queryLine.text.getBytes();

                    if (queryBytes.length != refBytes.length) {
                        throw new RuntimeException("Query and ref sequence unequal length");
                    }

                    String cigarString = "";
                    int nm = 0;
                    boolean indel = false;
                    for (int i = 0; i < queryBytes.length; i++) {
                        byte q = queryBytes[i];
                        byte ref = refBytes[i];

                        if (q == '-') {
                            if (ref == '-') {
                                //ignore, caused by insertion in another alignment;
                            } else {
                                cigarString += "D";
                                if (!indel) {
                                    indel = true;
                                    nm++;
                                }
                            }
                        } else {
                            if (ref == '-') {
                                cigarString += "I";
                                if (!indel) {
                                    indel = true;
                                    nm++;
                                }
                            } else {
                                indel = false;
                                cigarString += "M";
                                if (queryBytes[i] != refBytes[i]) {
                                    nm++;
                                }
                            }
                        }
                    }

                    cigarString = collapseCigar(cigarString);

                    String readName = queryLine.src;

                    String barcode = null;
                    if (readName.contains("Barcode")) {

                        int bcIdx = readName.indexOf("Barcode") + 8;
                        int endIdx = readName.indexOf(':', bcIdx);
                        barcode = readName.substring(bcIdx, endIdx);
                    }


                    // Output SAM record
                    int start = referenceLine.start + 1;
                    int mapq = 30;
                    String rnext = "*";
                    int pnext = 0;
                    int tlen = 0;
                    String seq = collapseSequence(queryLine.text);
                    String qual = "*";

                    out.print(readName + "\t" + flags + "\t" + chr + "\t" + start + "\t" + mapq + "\t" + cigarString + "\t" +
                            rnext + "\t" + pnext + "\t" + tlen + "\t" + seq + "\t" + qual + "\tNM:i:" + nm);
                    if (barcode != null) {
                        out.print("\tBC:Z:" + barcode);
                    }
                    out.println();
                }
            }
        }

        // Output SAN header
    }

    private static void outputHeader(Map<String, Integer> sequenceDictionary, PrintWriter out) {

        List<String> chrNames = new ArrayList<String>(sequenceDictionary.keySet());
        Collections.sort(chrNames, ChromosomeNameComparator.get());

        out.println("@HD\tVN:1.5");
        for (String chr : chrNames) {
            out.println("@SQ\tSN:" + chr + "\tLN:" + sequenceDictionary.get(chr));
        }

        out.println("@PG\tPN:MAFtoSAM\tID:MAFtoSAM");
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


    private static String collapseSequence(String text) {
        return text.replaceAll("-", "");
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

