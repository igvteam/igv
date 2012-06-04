package org.broad.igv.tools.converters;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.feature.Strand;
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

        System.out.println("BamToBed.  Pairs=" + properPairs);
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

            System.out.println("Max insert size = " + maxInsertSize);
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
