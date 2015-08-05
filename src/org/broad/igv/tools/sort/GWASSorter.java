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

package org.broad.igv.tools.sort;

import org.apache.log4j.Logger;
import org.broad.igv.gwas.GWASParser;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 4/12/12
 */
public class GWASSorter extends Sorter {

    static Logger log = Logger.getLogger(GWASSorter.class);

    List<String> headerLines = new ArrayList<String>();
    private GWASParser.GWASColumns columns;


    public GWASSorter(File inputFile, File outputFile) {
        super(inputFile, outputFile);
        this.columns = new GWASParser.GWASColumns();
        findColumns(inputFile);
    }

    private void findColumns(File inputFile) {
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(inputFile));
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                headerLines.add(nextLine);
                if (!nextLine.startsWith("#")) {
                    this.columns.parseHeader(nextLine);
                    break;
                }
            }
            if (this.columns.chrCol < 0) {
                throw new RuntimeException("Could not find chromosome column");
            }
            if (this.columns.locationCol < 0) {
                throw new RuntimeException("Could not find start column");
            }
        } catch (IOException e) {
            log.error("Error reading GWAS file", e);
            throw new RuntimeException("Error reading GWAS file" + e.getMessage(), e);
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {

            }
        }
    }


    @Override
    Parser getParser() {
        return new Parser(this.columns.chrCol, this.columns.locationCol, true);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException {

        String nextLine;
        while((nextLine = reader.readLine()) != null) {
            writer.println(nextLine);
            if (!nextLine.startsWith("#")) {
                break;
            }
        }
        return null;

    }
}
