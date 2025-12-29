package org.igv.tools.converters;

import org.igv.util.ParsingUtils;

import java.io.*;


/**
 * Created by jrobinson on 5/5/16.
 * <p>
 * Parses a MAF file for Alignment records, as opposed to Multiple Alignments.  For LAST aligner output.
 */

public class MFastaToSAM {

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
        String inputPath = args[5];
        String outputPath = args[6];
        convert(inputPath, outputPath);
    }


    public static void convert(String path, String outputPath) throws IOException {

        File outputFile = new File(outputPath);
        BufferedReader reader = null;
        PrintWriter out = null;

        try {
            reader = ParsingUtils.openBufferedReader(path);
            out = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            String chr = "MSA_SARS2_20200329_consensus";
            int chrLen = 30005;
            out.println("@HD\tVN:1.5");
            out.println("@SQ\tSN:" + chr + "\tLN:" + chrLen);
            out.println("@PG\tPN:MFastaToSAM\tID:MFastaToSAM");

            String readname = null;
            String seq = null;
            String cigarString = null;
            String line;
            int m = 0;
            int d = 0;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (readname != null) {
                        outputSamRecord(readname, chr, cigarString, seq, out);
                    }
                    readname = line.substring(1);
                    cigarString = "";
                    seq = "";
                    m = 0;
                    d = 0;
                } else {
                    char[] chars = line.toCharArray();
                    boolean lastCharDash = chars[0] == '-';
                    for (int i = 0; i < chars.length; i++) {
                        char c = chars[i];
                        if (c == '-') {
                            if (!lastCharDash) {
                                cigarString += String.valueOf(m) + 'M';
                                m = 0;
                            }
                            d++;
                            lastCharDash = true;
                        } else {
                            if (lastCharDash) {
                                cigarString += String.valueOf(d) + 'D';
                                d = 0;
                            }
                            m++;
                            seq += c;
                            lastCharDash = false;
                        }
                    }
                }
            }
            out.flush();

        } finally {
            out.close();
            reader.close();
        }
    }

    private static void outputSamRecord(String readName, String chr, String cigarString, String seq, PrintWriter out) {
        // Output SAM record
        out.print(readName + "\t" + 0 + "\t" + chr + "\t" + 1 + "\t" + 32 + "\t" + cigarString + "\t" +
                "*" + "\t" + 0 + "\t" + 0 + "\t" + seq + "\t*");

        out.println();
    }


}

