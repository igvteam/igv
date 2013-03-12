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

package org.broad.igv.feature.tribble;

import org.broad.igv.sam.reader.SamUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.ParsingUtils;

/**
 * NOTE: I don't think this class is used anymore, I think we just use picard for
 * all things SAM related.
 * Author: jrobinso
 * Date: Jul 18, 2010
 * Time: 3:48:18 PM
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
