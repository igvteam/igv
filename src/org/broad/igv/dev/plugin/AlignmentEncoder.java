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

package org.broad.igv.dev.plugin;

import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SAMWriter;
import org.broad.igv.sam.SamAlignment;

/**
 * Encode an Alignment into SAM format
 * User: jacob
 * Date: 2012-Sep-27
 */
public class AlignmentEncoder implements FeatureEncoder<Alignment> {

    //private boolean headerWritten = false;

    @Override
    public String encode(Alignment feature) {
        if (feature instanceof SamAlignment) {
            SamAlignment alignment = (SamAlignment) feature;
            String out = "";
            //TODO This is a hack, but in theory should work. However, we may be better off just
            //not having a header at all
//            if(!headerWritten){
//                out = alignment.getRecord().getHeader().getTextHeader() + "\n";
//                headerWritten = true;
//            }
            out += alignment.getRecord().getSAMString();
            return out;
        }
        return SAMWriter.getSAMString(feature);
    }

    @Override
    public int getNumCols(String line) {
        return line.split("\\t").length;
    }

    @Override
    public String getHeader() {
        return null;
    }
}
