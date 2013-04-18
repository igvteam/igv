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

import org.broad.igv.feature.tribble.MUTCodec;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author jrobinso
 *         Date: 4/8/13
 *         Time: 4:23 PM
 */
public class MUTSorter extends Sorter {

    int chrCol;
    int startCol;


    public MUTSorter(File inputFile, File outputFile) {

        super(inputFile, outputFile);
        MUTCodec codec = new MUTCodec(inputFile.getAbsolutePath(), null);
        chrCol = codec.getChrColumn();
        startCol = codec.getStartColumn();
    }


    @Override
    Parser getParser() throws IOException {

        return new Parser(chrCol, startCol);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException {

        String nextLine;
        // TODO -- check "browser" line syntax,  is it a multi-line directive?
        while ((nextLine = reader.readLine()).startsWith("#")) {
            writer.println(nextLine);   // Comments, directives
        }
        writer.println(nextLine);  // Header

        // Contract is to return the first data line!
        nextLine = reader.readLine();

        return nextLine;


    }
}
