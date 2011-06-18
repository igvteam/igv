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

package org.broad.igv.util;

import java.io.*;


/**
 * @author jrobinso
 */
public class ConvertPlasmodium {

    public static void main(String[] args) throws IOException {
        String inputDir = "/Users/jrobinso/plasmodium/Moran";

        for (File f : (new File(inputDir)).listFiles()) {

            if (f.getName().endsWith("wig")) {
                String outputName = f.getName() + ".txt";

                BufferedReader r = new BufferedReader(new FileReader(f));
                PrintWriter pw = new PrintWriter(new FileWriter(new File(f.getParentFile(), outputName)));

                String nextLine = null;
                while ((nextLine = r.readLine()) != null) {
                    pw.println(nextLine.replace("chrom=chr", "chrom=MAL"));
                }

                pw.close();
                r.close();
            }
        }

    }
}

