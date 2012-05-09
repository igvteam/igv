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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam.reader;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.DotAlignedAlignment;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.IOException;

/**
 * @author jrobinso
 */
public class DotAlignedParser implements AlignmentParser {

    private static Logger log = Logger.getLogger(GeraldParser.class);

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
