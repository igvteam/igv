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
