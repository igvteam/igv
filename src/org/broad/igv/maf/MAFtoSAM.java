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
 * <p/>
 * Parses a MAF file for Alignment records, as opposed to Multiple Alignments.  For LAST aligner output.
 */

public class MAFtoSAM {


    public static void main(String[] args) throws IOException {
        String inputPath = args[0];
        String outputPath = args.length > 1 ? args[1] : null;
        convert(inputPath, outputPath, true, false);
    }


    public static void convert(String path, String outputPath, boolean includeSequence, boolean groupByReadName) throws IOException {

        String samOutput = outputPath.endsWith(".sam") ? outputPath : outputPath + ".sam";
        String unsortedOutput = outputPath + ".unsorted.sam";
        String sortedOutput = outputPath + "sorted.sam";
        int readCounter = 1;

        BufferedReader reader = ParsingUtils.openBufferedReader(path);
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(unsortedOutput)));
        Map<String, Integer> sequenceDictionary = new LinkedHashMap<String, Integer>();
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("a ")) {
                // Parse alignments until blank line
                parseBlock(reader, out, sequenceDictionary, includeSequence);
            }
        }
        out.flush();
        out.close();
        reader.close();

        if (groupByReadName) {
            unsortedOutput = combineReads(unsortedOutput);
        }

        // Now sort sam records
        Sorter sorter = SorterFactory.getSorter(new File(unsortedOutput), new File(sortedOutput));
        sorter.run();

        // Finally insert sam header
        out = new PrintWriter(new BufferedWriter(new FileWriter(samOutput)));
        outputHeader(sequenceDictionary, out);

        reader = ParsingUtils.openBufferedReader(sortedOutput);
        while ((line = reader.readLine()) != null) {
            out.println(line);
        }
        out.flush();
        out.close();
        reader.close();

        (new File(unsortedOutput)).deleteOnExit();
        (new File(sortedOutput)).deleteOnExit();

    }

    private static String combineReads(String unsortedOutput) throws IOException {

        String groupedOutput = unsortedOutput + ".grouped.sam";
        String combinedOutput = unsortedOutput + ".combined.sam";

        Sorter sorter = SorterFactory.getSorter(new File(unsortedOutput), new File(groupedOutput));
        sorter.setComparator(SAMSorter.ReadNameComparator);
        sorter.run();

        // Get Sam file reader, output header, then iterate through sam records combining those with like read names

        BufferedReader reader = new BufferedReader(new FileReader(groupedOutput));
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(combinedOutput)));
        String line;
        String[] lastRecord = null;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("@")) {
                out.println(line);
            } else {
                String[] record = Globals.tabPattern.split(line);

                if (lastRecord != null && lastRecord[0].equals(record[0])) {

                    Cigar lastCigar = TextCigarCodec.decode(lastRecord[5]);
                    int lastEnd = Integer.parseInt(lastRecord[3]) + lastCigar.getReferenceLength()-1;
                    int gap = Integer.parseInt(record[3]) - lastEnd - 1;

                    if(gap < 1) {
                     // Don't try to combine overlapping records
                        printRecord(out, lastRecord);
                        lastRecord = record;
                    }
                    else {

                        String newCigar = lastRecord[5] + (gap + "D") + record[5];
                        lastRecord[5] = newCigar;
                        if (lastRecord[9].equals("*") || record[9].equals("*")) {
                            lastRecord[9] = "*";
                        } else {
                            lastRecord[9] = lastRecord[9] + record[9];
                        }
                    }

                } else {
                    if(lastRecord != null) printRecord(out, lastRecord);
                    lastRecord = record;
                }
            }
        }

        // Last one
        if(lastRecord != null) printRecord(out, lastRecord);

        out.flush();
        out.close();

       // (new File(groupedOutput)).deleteOnExit();
        (new File(combinedOutput)).deleteOnExit();

        return combinedOutput;

    }

    private static void printRecord(PrintWriter out, String [] record) {
        out.print(record[0]);
        for(int i=1; i<record.length; i++) {
            out.print("\t" + record[i]);
        }
        out.println();
    }


    private static SAMRecord combineRecords(SAMRecord lastRecord, SAMRecord record) {
        return null;
    }

    private static void parseBlock(BufferedReader reader, PrintWriter out, Map<String, Integer> sequenceDictionary,
                                   boolean includeSequence) throws IOException {

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

                    // Build cigar string one byte at a time
                    byte[] queryBytes = queryLine.text.getBytes();

                    if (queryBytes.length != refBytes.length)
                        throw new RuntimeException("Query and ref sequence unequal length");

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
                    String qname = queryLine.src;
                    int flag = 0;

                    int start = referenceLine.start + 1;
                    int mapq = 30;
                    String rnext = "*";
                    int pnext = 0;
                    int tlen = 0;
                    String seq = includeSequence ? collapseSequence(queryLine.text) : "*";
                    String qual = "*";

                    out.println(qname + "\t" + flag + "\t" + chr + "\t" + start + "\t" + mapq + "\t" + cigarString + "\t" +
                            rnext + "\t" + pnext + "\t" + tlen + "\t" + seq + "\t" + qual);
                }
            }
        }

        // Output SAN header
    }

    private static void outputHeader(Map<String, Integer> sequenceDictionary, PrintWriter out) {

        List<String> chrNames = new ArrayList<String>(sequenceDictionary.keySet());
        Collections.sort(chrNames, ChromosomeNameComparator.get());

        out.println("@HD\tVN:1.5\t SO:coordinate");
        for (String chr : chrNames) {
            out.println("@SQ\tSN:" + chr + "\tLN:" + sequenceDictionary.get(chr));
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

