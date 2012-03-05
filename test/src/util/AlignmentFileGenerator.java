package util;

import org.junit.Ignore;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Generates alignment files (.aligned format) for testing
 *
 * @author Jim Robinson
 * @date 3/2/12
 */
@Ignore
public class AlignmentFileGenerator {

    public static void main(String[] args) throws IOException {

        String outputFile = "test.aligned";
        String chr = "chr1";
        int averageDepth = 1000;
        int readSize = 100;
        int interval = 1000;

        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new FileWriter(outputFile));

            int averageReadsPerPosition = averageDepth / readSize;

            int startPosition = 100; // arbitrary
            for (int pos = startPosition; pos < startPosition + interval; pos++) {

                int readCount = (int) (averageReadsPerPosition * 2 * Math.random());
                while (readCount-- > 0) {
                    pw.println(chr + "\t" + pos + "\t" + (pos + readSize) + "\t+");
                }

            }

        } finally {
            if (pw != null) pw.close();
        }

    }
}
