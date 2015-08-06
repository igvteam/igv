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
