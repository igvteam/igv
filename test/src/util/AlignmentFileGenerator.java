package util;

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

    public static void main(String[] args) throws IOException {

        String outputFile = "test.aligned";
        String chr = "chr1";
        int averageDepth = 1000;
        int avgReadSize = 50;
        int readSigma = 5;
        int interval = 200;
        boolean append = false;

        int startPosition = 800; // arbitrary

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
