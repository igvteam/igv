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

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;

import java.io.*;

/**
 * Converts a bam -> a bed file by writing each record as a bed feature.
 * <p/>
 * If the "properPair" option is true paired alignments marked "proper" are output as a single record, with the
 * name field containing the insert size.
 *
 * @author Jim Robinson
 * @date 6/4/12
 */
public class BamToBed {

    public static void convert(File inputBam, File outputBed, boolean properPairs) throws IOException {

        AlignmentReader reader = null;
        CloseableIterator<Alignment> iter = null;
        PrintWriter bedWriter = null;

        int maxInsertSize = 0;

        try {
            bedWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputBed)));
            reader = AlignmentReaderFactory.getReader(inputBam.getAbsolutePath(), false);
            iter = reader.iterator();

            while (iter.hasNext()) {
                Alignment a = iter.next();
                if(passFilter(a, properPairs)) {
                    int start = a.getAlignmentStart();
                    final int insertSize = Math.abs(a.getInferredInsertSize());
                    int end = properPairs ? (start + insertSize) : a.getAlignmentEnd();
                    String name = a.getReadName();
                    String strand = properPairs ? "." : (a.isNegativeStrand() ? "-" : "+");
                    bedWriter.print(a.getChr() + "\t" + start + "\t" + end + "\t" + name + "\t" + strand);
                    if(properPairs) {
                        bedWriter.println("\t" + insertSize);
                    }
                    else {
                        bedWriter.println();
                    }

                    maxInsertSize = insertSize > maxInsertSize ? insertSize : maxInsertSize;
                }

            }
        } finally {
            if(bedWriter != null) bedWriter.close();
            if(iter != null) iter.close();
            if(reader != null) reader.close();
        }
    }



    private static boolean passFilter(Alignment alignment, boolean properPairs) {

         // For paired coverage, see if the alignment is properly paired, and if it is the "leftmost" alignment
        // (to prevent double-counting the pair).
        if(properPairs) {
            ReadMate mate = alignment.getMate();
            if(!alignment.isProperPair() || alignment.getMate() == null || alignment.getStart() > mate.getStart()) {
                return false;
            }
        }

        return alignment.isMapped()  && !alignment.isDuplicate() && !alignment.isVendorFailedRead();
    }

}
