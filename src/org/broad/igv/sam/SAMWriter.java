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

package org.broad.igv.sam;

import net.sf.samtools.*;

import java.io.File;
import java.io.OutputStream;

/**
 * Write SAM/BAM Alignments to a file or stream
 * <p/>
 * User: jacob
 * Date: 2012/05/04
 */
public class SAMWriter {

    private SAMFileHeader header;

    public SAMWriter(SAMFileHeader header) {
        this.header = header;
    }

    public void writeToFile(File outFile, Iterable<SamAlignment> alignments) {
        SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, outFile);
        writeAlignments(writer, alignments);
    }

    public void writeToStream(OutputStream stream, Iterable<SamAlignment> alignments, boolean bam) {

        SAMFileWriterImpl writer;
        if (bam) {
            //BAMFileWriter can't take null argument for File.
            //Have sent in a patch. May 4 2012
            writer = new BAMFileWriter(stream, new File(""));
        } else {
            writer = new SAMTextWriter(stream);
        }

        writer.setHeader(header);
        writeAlignments(writer, alignments);
    }

    private void writeAlignments(SAMFileWriter writer, Iterable<SamAlignment> alignments) {
        for (SamAlignment alignment : alignments) {
            writer.addAlignment(alignment.getRecord());
        }
        writer.close();
    }
}
