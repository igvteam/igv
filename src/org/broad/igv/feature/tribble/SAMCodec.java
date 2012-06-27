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

package org.broad.igv.feature.tribble;

import org.broad.igv.sam.reader.SamUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.ParsingUtils;

/**
 * Author: jrobinso
 * Date: Jul 18, 2010
 * Time: 3:48:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class SAMCodec extends AsciiFeatureCodec {

    final static int FLAG_COL = 1;
    final static int READ_UNMAPPED_FLAG = 0x4;

    // From SAM specification

    private static final int RNAME_COL = 2;
    private static final int POS_COL = 3;
    private static final int CIGAR_COL = 5;
    private static final int NUM_REQUIRED_FIELDS = 11;

    // Declare static for effeciency -- note THIS IS NOT THREAD SAFE
    private static String [] fields = new String[NUM_REQUIRED_FIELDS];

    protected SAMCodec() {
        super(SAMCodec.class);
    }

    public Feature decodeLoc(String line) {
        return decode(line);  //To change body of implemented methods use File | Settings | File Templates.
    }

    public synchronized Feature decode(String line) {

        if(line.startsWith("@")) {
            return null;
        }

        final int numFields = ParsingUtils.split(line, fields, '\t');
        if (numFields < NUM_REQUIRED_FIELDS ) {
            return null;
        }

        String chr = fields[RNAME_COL];
        int alignmentStart = Integer.parseInt(fields[POS_COL].trim());
        String cigarString = fields[CIGAR_COL];
        int len = SamUtils.getPaddedReferenceLength(cigarString);
        int alignmentEnd = alignmentStart + len - 1;

        return new Locus(chr, alignmentStart, alignmentEnd);

    }

    public Object readHeader(LineReader reader) {
        return null;
    }

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     * <p/>
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     * <p/>
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(String path) {
        return path.toLowerCase().endsWith(".sam");
    }

    boolean isMapped(String[] fields) {
        int flags = Integer.parseInt(fields[FLAG_COL]);
        return (flags & READ_UNMAPPED_FLAG) == 0;
    }
}
