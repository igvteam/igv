
package org.igv.maf;

import org.igv.util.ParsingUtils;

import java.io.*;
import java.util.Properties;

/**
 * @author jrobinso
 */
public class MAFUtils {

    public static void main(String[] args) throws IOException {
        String maf = "/Users/jrobinso/projects/maf/danRer7.gasAcu1.net.maf";
        createTestFile(maf, 100, 100);

    }

    /**
     * Create a test file by keeping a sampling of blocks from the input file
     *
     * @param mafFile
     * @param sampling
     * @param maxBlocks
     */
    public static void createTestFile(String mafFile, int sampling, int maxBlocks) throws IOException {

        BufferedReader reader = new BufferedReader(new FileReader(mafFile));
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("Test.maf")));

        int a = 0;
        int n = 0;

        String line;
        while ((line = reader.readLine()) != null && n < maxBlocks) {
            if (line.startsWith("#")) {
                pw.println(line);
            }


            String[] tokens = ParsingUtils.WHITESPACE_PATTERN.split(line);
            if (tokens[0].equals("a")) {
                // Peek
                    String tmp = reader.readLine();
                tokens = ParsingUtils.WHITESPACE_PATTERN.split(tmp);


                if (a % sampling == 0 && tokens[1].contains("chr")) {
                    pw.println();
                    pw.println(line);
                    pw.println(tmp);
                    while ((line = reader.readLine()) != null) {
                        tokens = ParsingUtils.WHITESPACE_PATTERN.split(line);
                        if (tokens[0].equals("a")) {
                            break;
                        } else {
                            pw.println(line);
                        }

                    }
                }

                a++;
            }
       }

        pw.close();
        reader.close();

    }
}
