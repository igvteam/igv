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

package org.broad.igv.util.converters;

import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.List;

/**
 * Converts a directory of "density" files (Tarjei et al) to bedgraph
 *
 * @author Jim Robinson
 * @date 9/13/11
 */
public class DensitiesToBedGraph {


    public static void main(String[] args) throws IOException {


        File inputDir = new File(args[0]);
        File outputDir = new File(args[1]);

        convertAll(inputDir, outputDir);

    }

    private static void convertAll(File inputDir, File outputDir) throws IOException {
        File[] files = inputDir.listFiles();
        for (File f : files) {
            if (f.getAbsolutePath().endsWith(".densities.txt.gz")) {
                String ofile = f.getName().replace(".densities.txt.gz", ".bedgraph");
                convert(f, new File(outputDir, ofile));
            }
        }
    }

    public static void convert(File ifile, File ofile) throws IOException {

        BufferedReader reader = null;
        PrintWriter pw = null;

        reader = ParsingUtils.openBufferedReader(ifile.getAbsolutePath());
        pw = new PrintWriter(new BufferedWriter(new FileWriter(ofile)));

        String nextLine;
        while ((nextLine = reader.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            float counts = Float.parseFloat(tokens[2]);
            if (counts > 0) {
                int start = Integer.parseInt(tokens[1]);
                pw.println(tokens[0] + "\t" + start + "\t" + (start + 25) + "\t" + counts);
            }
        }

        reader.close();
        pw.close();
    }
}
