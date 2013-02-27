/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
