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

package org.broad.igv.tools.sort;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ParsingUtils;

import java.io.IOException;

/**
 * @author: nazaire
 */
public class Parser {

    private int chrCol;
    private int startCol;
    private String[] fields;

    public Parser(int chrCol, int startCol) {
        this.chrCol = chrCol;
        this.startCol = startCol;
        fields = new String[this.startCol + 1];
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
        }

        return createRecord(nextLine);
    }

    public SortableRecord createRecord(String nextLine) {
        int nTokens = ParsingUtils.split(nextLine, fields, '\t');
        // TODO -- what to do if nTokens < startCol?

        String chr = fields[chrCol];

        int start ;
        try {
            start = Integer.parseInt(fields[startCol].trim());
        }
        catch(NumberFormatException e) {
            start = Integer.MAX_VALUE;
        }
        String text = nextLine;

        return new SortableRecord(chr, start, text);
    }
}
