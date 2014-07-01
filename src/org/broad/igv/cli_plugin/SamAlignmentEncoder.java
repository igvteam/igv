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

package org.broad.igv.cli_plugin;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriterImpl;
import net.sf.samtools.SAMTextWriter;
import org.broad.igv.sam.PicardAlignment;

import java.io.OutputStream;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Encode an Alignment into SAM format
 * @author jacob
 * @since 2012-Sep-27
 */
public class SamAlignmentEncoder implements FeatureEncoder<PicardAlignment> {

    private boolean headerSet = false;

    public Map<String, Object> encodeAll(OutputStream stream, Iterator<? extends PicardAlignment> alignments) {
        SAMFileWriterImpl writer = new SAMTextWriter(stream);
        while (alignments.hasNext()) {
            PicardAlignment samAl = alignments.next();
            if (!headerSet) {
                writer.setSortOrder(SAMFileHeader.SortOrder.unsorted, true);
                writer.setHeader(samAl.getRecord().getHeader());
                headerSet = true;
            }
            writer.addAlignment(samAl.getRecord());
        }
        writer.close();
        return null;
    }



    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap, Argument argument) {
        //pass
    }
}
