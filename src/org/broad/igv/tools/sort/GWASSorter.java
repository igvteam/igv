/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.tools.sort;

import org.apache.log4j.Logger;
import org.broad.igv.gwas.GWASParser;
import org.broad.tribble.readers.AsciiLineReader;

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
