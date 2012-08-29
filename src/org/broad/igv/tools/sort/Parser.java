/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.tools.sort;

import org.broad.igv.Globals;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.IOException;

/**
 * @author: nazaire
 */
public class Parser {

    private int chrCol;
    private int startCol;
    boolean splitOnWhiteSpace;
    private String commentPrefix;

    public Parser(int chrCol, int startCol) {
        this(chrCol, startCol, false, "#");
    }

    public Parser(int chrCol, int startCol, boolean splitOnWhiteSpace) {
        this(chrCol, startCol, splitOnWhiteSpace, "#");
    }

    public Parser(int chrCol, int startCol, boolean splitOnWhitespace, String commentPrefix) {
        this.chrCol = chrCol;
        this.startCol = startCol;
        this.splitOnWhiteSpace = splitOnWhitespace;
        this.commentPrefix = commentPrefix;
    }


    public SortableRecord readNextRecord(AsciiLineReader reader) {
        String nextLine = null;
        try {
            nextLine = reader.readLine();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            return null;
        }
        if (nextLine == null) {
            return null;
        } else if (nextLine.startsWith(commentPrefix)) {
            return readNextRecord(reader);
        }

        return createRecord(nextLine);
    }

    public SortableRecord createRecord(String nextLine) {
        String[] fields = splitOnWhiteSpace ?
                Globals.singleTabMultiSpacePattern.split(nextLine) :
                Globals.tabPattern.split(nextLine, -1);

        String chr = fields[chrCol];

        int start;
        try {
            start = Integer.parseInt(fields[startCol].trim());
        } catch (NumberFormatException e) {
            start = Integer.MAX_VALUE;
        }
        String text = nextLine;

        return new SortableRecord(chr, start, text);
    }
}
