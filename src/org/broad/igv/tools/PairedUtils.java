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

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.feature.Strand;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by jrobinso on 10/6/15.
 */
public class PairedUtils {

    public static void main(String [] args) {
        extractInteractions("/Users/jrobinso/Downloads/wgEncodeGisDnaPetK562F10kAln.bam", "test.interactions", 1000);
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

                if (passFilter(alignment)) {

                    String chr1 = alignment.getChr();
                    int bin1 = alignment.getAlignmentStart() / binSize;
                    String chr2 = alignment.getMate().getChr();
                    int bin2 = alignment.getMate().getStart() / binSize;

                    String cell =chr1 + "_" +  bin1 + ":" + chr2 + "_" + bin2;

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
            while(iter2.hasNext()) {
                if(iter2.next().getValue() < 10) iter2.remove();
            }

            // Output counts

            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            for(Map.Entry<String, Integer> entry : entrySet) {

                String cell =  entry.getKey();
                String [] tokens = cell.split(":");
                String [] t1 = tokens[0].split("_");
                String [] t2 = tokens[1].split("_");

                int bin1 = Integer.parseInt(t1[1]);
                int bin2 = Integer.parseInt(t2[1]);
                int gPos1 = (int) (bin1 + 0.5) * binSize;
                int gPos2 = (int) (bin2 + 0.5) * binSize;
                pw.println(cell + "\t" + t1[0] + "\t" + gPos1 + "\t" + t2[0] + "\t" + gPos2 + "\t" + (entry.getValue()/2));
            }


        } catch (Exception e) {
            e.printStackTrace();
        } finally {

            if(pw != null) {
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

    private static boolean passFilter(Alignment alignment) {
        return alignment.isPaired() && alignment.getMate().isMapped() &&
                (Math.abs(alignment.getInferredInsertSize()) > 100000 || !alignment.getChr().equals(alignment.getMate().getChr()));
    }

}
