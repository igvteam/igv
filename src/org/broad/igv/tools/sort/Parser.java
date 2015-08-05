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

import org.broad.igv.Globals;
import htsjdk.tribble.readers.AsciiLineReader;

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
            e.printStackTrace();
            return null;
        }
        if (nextLine == null) {
            return null;
        } else if (nextLine.trim().length() == 0 || nextLine.startsWith(commentPrefix)) {
            return readNextRecord(reader);
        }

        try {
            return createRecord(nextLine);
        } catch (ArrayIndexOutOfBoundsException e) {
            e.printStackTrace();
            System.out.println(nextLine);
            throw e;
        }
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
