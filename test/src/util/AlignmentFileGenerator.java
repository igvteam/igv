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

        String outputFile = "wide_spread.aligned";
        String chr = "chr1";
        int averageDepth = 2000;
        int avgReadSize = 100;
        int readSigma = 20;
        int interval = 1000;

        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new FileWriter(outputFile));

            int averageReadsPerPosition = averageDepth / avgReadSize;

            int startPosition = 100; // arbitrary
            for (int pos = startPosition; pos < startPosition + interval; pos++) {
                
                int readSize = avgReadSize + (int) (readSigma * RAND.nextGaussian());
                int readCount = RAND.nextInt(2*averageReadsPerPosition);
                while (readCount-- > 0) {
                    pw.println(chr + "\t" + pos + "\t" + (pos + readSize) + "\t+");
                }

            }

        } finally {
            if (pw != null) pw.close();
        }

    }
}
