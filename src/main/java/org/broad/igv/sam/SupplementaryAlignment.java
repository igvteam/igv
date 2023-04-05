package org.broad.igv.sam;

import htsjdk.samtools.util.Locatable;
import org.broad.igv.Globals;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class SupplementaryAlignment implements Locatable {

    public String chr;
    public int start;
    public char strand;
    public int mapQ;
    public int numMismatches;
    public int lenOnRef;


    public SupplementaryAlignment(String rec) {
        String[] tokens = Globals.commaPattern.split(rec);
        chr = tokens[0];
        start = Integer.parseInt(tokens[1]);
        strand = tokens[2].charAt(0);
        mapQ = Integer.parseInt(tokens[4]);
        numMismatches = Integer.parseInt(tokens[5]);
        lenOnRef = computeLengthOnReference(tokens[3]);
    }

    public String printString() {
        // chr6:43,143,415-43,149,942 (-) @ MAPQ 60 NM 763
        return chr + ":" + Globals.DECIMAL_FORMAT.format(start) + "-" + Globals.DECIMAL_FORMAT.format(start + lenOnRef)
                + " (" + strand + ") = " + Globals.DECIMAL_FORMAT.format(lenOnRef) + "bp  @MAPQ " + mapQ + " NM" + numMismatches;
    }

    public static List<SupplementaryAlignment> parseFromSATag(String saTag){
        return Arrays.stream(Globals.semicolonPattern.split(saTag))
                .map(SupplementaryAlignment::new)
                .collect(Collectors.toList());
    }

    int computeLengthOnReference(String cigarString) {

        int len = 0;
        StringBuffer buf = new StringBuffer();

        for (char c : cigarString.toCharArray()) {

            if (c > 47 && c < 58) {
                buf.append(c);
            } else {
                switch (c) {
                    case 'N':
                    case 'D':
                    case 'M':
                    case '=':
                    case 'X':
                        len += Integer.parseInt(buf.toString());
                }
                buf.setLength(0);
            }

        }
        return len;
    }

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return start + lenOnRef;
    }
}
