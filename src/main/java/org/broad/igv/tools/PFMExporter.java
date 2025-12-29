package org.broad.igv.tools;

import org.broad.igv.sam.AlignmentCounts;
import org.broad.igv.util.StringUtils;

import java.util.ArrayList;

/**
 * Export the consensus of alignment counts. Format is similar to
 * http://lgsun.grc.nia.nih.gov/CisFinder/cisfinder-help.html#format
 * <p/>
 * We have 2 header lines at the top, followed by a matrix with the count
 * of each base in each column, and the position along rows
 * Our format is similar but the header is missing many optional fields
 * > [locus]    [consensus]
 * > A  C   G   T   N
 * #a   #c  #g  #t #n
 * ....
 * ...
 *
 * @author jacob
 * @date 2013-Jul-01
 */
public class PFMExporter {

    private static final char[] nucleotides = {'a', 'c', 'g', 't', 'n'};
    private static final int numNuc = nucleotides.length;

    public static String header1;

    private static String delim = "\t";

    static {
        String[] capNuc = new String[numNuc];
        for (int n = 0; n < capNuc.length; n++) {
            capNuc[n] = String.valueOf(nucleotides[n]).toUpperCase();
        }
        header1 = StringUtils.join(capNuc, delim);
    }

    /**
     * @param counts
     * @param chr
     * @param start  0-based
     * @param end    0-based, end-exclusive
     * @return
     */
    public static String createPFMText(AlignmentCounts counts, String chr, int start, int end) {

        AbstractConsensusCalculator consCalc = new CavenerConsensusCalculator();
        ArrayList<String> allTextList = new ArrayList<String>(2);


        char[] consensus = new char[end - start + 1];
        for (int loc = start; loc < end; loc++) {
            consensus[loc - start] = consCalc.calculateConsensusBase(counts, loc);

            String[] curLocCounts = new String[numNuc];
            for (int ii = 0; ii < numNuc; ii++) {
                char c = nucleotides[ii];
                int negCount = counts.getNegCount(loc, (byte) c);
                int posCount = counts.getPosCount(loc, (byte) c);
                curLocCounts[ii] = String.format("%d", negCount + posCount);
            }
            String curLine = StringUtils.join(curLocCounts, delim);
            allTextList.add(curLine);
        }

        String topLine = String.format("%s:%d-%d%s%s", chr, start + 1, end, delim, String.valueOf(consensus).toUpperCase());

        //Add in reverse order because adding in beginning
        allTextList.add(0, header1);
        allTextList.add(0, topLine);

        return StringUtils.join(allTextList, "\n");

    }
}
