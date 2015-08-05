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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.maf;

import org.broad.igv.util.ParsingUtils;

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
