/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package util;

import org.broad.igv.util.TestUtils;
import org.junit.Ignore;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

/**
 * Generates alignment files (.aligned format) for testing
 *
 * @author Jim Robinson
 * @date 3/2/12
 */
@Ignore
public class AlignmentFileGenerator {

    static Random RAND = new Random();
    private static final int READ_PAIRED_FLAG = 0x1;
    private static final int PROPER_PAIR_FLAG = 0x2;
    private static final int READ_UNMAPPED_FLAG = 0x4;
    private static final int MATE_UNMAPPED_FLAG = 0x8;
    private static final int READ_STRAND_FLAG = 0x10;
    private static final int MATE_STRAND_FLAG = 0x20;
    private static final int FIRST_OF_PAIR_FLAG = 0x40;
    private static final int SECOND_OF_PAIR_FLAG = 0x80;
    private static final int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    private static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    private static final int DUPLICATE_READ_FLAG = 0x400;

    public static void main(String[] args) {

        int t = 99;
        System.out.println("READ_PAIRED_FLAG = " + (t & READ_PAIRED_FLAG));

         t |= DUPLICATE_READ_FLAG;
        System.out.println(t);
    }

    public static void generateTestFile() throws IOException {

        String outputFile = TestUtils.DATA_DIR + "aligned/short_spread.aligned";
        String chr = "chr1";
        int averageDepth = 1000;
        int avgReadSize = 100;
        int readSigma = 5;
        int interval = 1000;
        boolean append = false;

        int startPosition = 100; // arbitrary

        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new FileWriter(outputFile, append));

            float averageReadsPerPosition = (float) averageDepth / avgReadSize;

            for (int pos = startPosition; pos < startPosition + interval; pos++) {

                int readSize = avgReadSize + (int) (readSigma * RAND.nextGaussian());
                float avgReadCount = 2 * averageReadsPerPosition * RAND.nextFloat();
                int readCount = (int) avgReadCount;
                readCount += avgReadCount - readCount < RAND.nextFloat() ? 1 : 0;
                while (readCount-- > 0) {
                    pw.println(chr + "\t" + pos + "\t" + (pos + readSize) + "\t+");
                }

            }

        } finally {
            if (pw != null) pw.close();
        }

    }
}
