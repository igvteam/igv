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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam.reader;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.DotAlignedAlignment;
import org.broad.igv.util.ParsingUtils;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;

/**
 * @author jrobinso
 */
public class DotAlignedParser implements AlignmentParser {

    private static Logger log = Logger.getLogger(DotAlignedParser.class);

    public static int CHROMOSOME_COLUMN = 0;
    private static int START_COLUMN = 1;
    private static int END_COLUMN = 2;
    private static int STRAND_COLUMN = 3;
    private static int NAME_COLUMN = -1;
    boolean bedFormat = false;

    public DotAlignedParser() {
        bedFormat = false;
    }

    public DotAlignedParser(boolean bedFormat) {
        this.bedFormat = bedFormat;
        if (bedFormat) {
            STRAND_COLUMN = 5;
            NAME_COLUMN = 3;
        }
    }

    public DotAlignedAlignment readNextRecord(AsciiLineReader reader) {
        String nextLine;
        try {
            while ((nextLine = reader.readLine()) != null) {
                DotAlignedAlignment alignment = createAlignment(nextLine);
                if (alignment != null) {
                    return alignment;
                }
            }
        } catch (IOException e) {
            log.error("Error reading line", e);
        }
        return null;
    }

    private DotAlignedAlignment createAlignment(String nextLine) {


        try {
            String[] fields = Globals.tabPattern.split(nextLine, -1);
            int nTokens = fields.length;

            if (nTokens <= END_COLUMN) {
                System.out.println("Skipping line: " + nextLine);
                return null;
            }

            String chr = fields[CHROMOSOME_COLUMN];
            int start = Integer.parseInt(fields[START_COLUMN]);
            int end = Integer.parseInt(fields[END_COLUMN]);


            if (bedFormat) {
                boolean isNegative = false;
                String name = "";
                if (nTokens > STRAND_COLUMN) {
                    isNegative = fields[STRAND_COLUMN].equals("-");
                }
                if (nTokens > NAME_COLUMN) {
                    name = fields[NAME_COLUMN];
                }
                return new DotAlignedAlignment(chr, start, end, isNegative, name);
            } else {

                boolean isNegative = fields[STRAND_COLUMN].equals("-");
                return new DotAlignedAlignment(chr, start, end, isNegative);
            }
        } catch (NumberFormatException e) {
            System.out.println("Skipping line: " + nextLine);
            return null;
        }
    }
}
