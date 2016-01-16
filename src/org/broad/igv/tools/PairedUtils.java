/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 UC San Diego
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

package org.broad.igv.tools;

import htsjdk.samtools.*;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.stream.IGVSeekableBufferedStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.*;
import java.util.*;

/**
 * Created by jrobinso on 10/6/15.
 */
public class PairedUtils {

    public static void main(String[] args) {

        extractInteractions(args[0], args[1], Integer.parseInt(args[2]));
    }

    public static void extractInteractions(String alignmentFile, String outputFile, int binSize) {


        AlignmentReader reader = null;
        PrintWriter pw = null;
        CloseableIterator<Alignment> iter = null;
        Map<String, Integer> counts = new HashMap<String, Integer>(10000);

        try {

            reader = AlignmentReaderFactory.getReader(alignmentFile, false);
            iter = reader.iterator();


            while (iter != null && iter.hasNext()) {

                Alignment alignment = iter.next();

                if (alignment.isPaired() && alignment.getMate().isMapped() &&
                        alignment.getMappingQuality() > 0 &&
                        interactionFilter(alignment)) {

                    String chr1 = alignment.getChr();
                    int bin1 = alignment.getAlignmentStart() / binSize;
                    String chr2 = alignment.getMate().getChr();
                    int bin2 = alignment.getMate().getStart() / binSize;

                    int o1 = getOrder(chr1);
                    int o2 = getOrder(chr2);
                    if (o1 > o2 || (o1 < 0 || o2 < 0))
                        continue;  // Only need to record one diagonal, they are symmetrical

                    String cell = chr1 + "_" + bin1 + ":" + chr2 + "_" + bin2;

                    Integer count = counts.get(cell);
                    if (count == null) {
                        counts.put(cell, 1);
                    } else {
                        count += 1;
                        counts.put(cell, count);
                    }
                }
            }

            iter.close();

            // Now filter cells with counts < 5
            Set<Map.Entry<String, Integer>> entrySet = counts.entrySet();
            Iterator<Map.Entry<String, Integer>> iter2 = entrySet.iterator();
            while (iter2.hasNext()) {
                if (iter2.next().getValue() < 5) iter2.remove();
            }

            // Output counts

            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            for (Map.Entry<String, Integer> entry : entrySet) {

                String cell = entry.getKey();
                String[] tokens = cell.split(":");
                String[] t1 = tokens[0].split("_");
                String[] t2 = tokens[1].split("_");

                int bin1 = Integer.parseInt(t1[1]);
                int bin2 = Integer.parseInt(t2[1]);
                int gPos1 = (int) (bin1 + 0.5) * binSize;
                int gPos2 = (int) (bin2 + 0.5) * binSize;
                pw.println(cell + "\t" + t1[0] + "\t" + gPos1 + "\t" + t2[0] + "\t" + gPos2 + "\t" + (entry.getValue() / 2));
            }


        } catch (Exception e) {
            e.printStackTrace();
        } finally {

            if (pw != null) {
                pw.close();
            }
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

        }


    }

    public static void extractFunnyPairs(String alignmentFile, String outputFile) throws IOException {


        htsjdk.samtools.SamReader reader = null;
        SAMFileWriter writer = null;

        try {

            final SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

            SeekableStream ss = new IGVSeekableBufferedStream(IGVSeekableStreamFactory.getInstance().getStreamFor(alignmentFile), 128000);
            SamInputResource resource = SamInputResource.of(ss);
            reader = factory.open(resource);


            SAMFileHeader header = reader.getFileHeader();
            boolean preSorted = true;
            SAMFileWriterFactory wfactory = new SAMFileWriterFactory();

            writer = wfactory.makeSAMOrBAMWriter(header, preSorted, new File(outputFile));


            SAMRecordIterator iter = reader.iterator();

            int count = 0;
            while (iter != null && iter.hasNext()) {

                SAMRecord record = iter.next();
                Alignment alignment = new PicardAlignment(record);
                if (funnyPairFilter(alignment)) {
                    writer.addAlignment(record);
                }

                count++;
                //if(count > 1000000) break;
            }

            iter.close();

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            writer.close();
            reader.close();
        }


    }

    private static int getOrder(String chr) {
        int order;
        try {
            order = Integer.parseInt(chr);
        } catch (NumberFormatException e) {
            if (chr.contains("X")) {
                order = 23;
            } else if (chr.contains("Y")) {
                order = 24;
            } else {
                order = -1;
            }
        }

        return order;
    }

    private static boolean interactionFilter(Alignment alignment) {
        return alignment.isPaired() && alignment.getMate().isMapped() &&
                (Math.abs(alignment.getInferredInsertSize()) > 100000 || !alignment.getChr().equals(alignment.getMate().getChr()));
    }

    private static boolean funnyPairFilter(Alignment alignment) {

        if (!(alignment.isPaired() && alignment.getMate().isMapped() && alignment.getMappingQuality() > 0)) return false;

        // Orientation hard-coded for FR
        if (!(frTypes.contains(alignment.getPairOrientation()))) return true;

        // Insert size, again hard-coded
        if (Math.abs(alignment.getInferredInsertSize()) > 1000) return true;

        // Mate chromosome
        if (alignment.getMate().isMapped() && !alignment.getChr().equals(alignment.getMate().getChr())) return true;

        // Mate unmapped
        // if (!alignment.getMate().isMapped()) return true;

        return false;

    }

    private static HashSet<String> frTypes = new HashSet<String>();
    static {
        frTypes.add("F1R2");
        frTypes.add("F2R1");
        frTypes.add("FR");
        frTypes.add("F R ");
    }


}
