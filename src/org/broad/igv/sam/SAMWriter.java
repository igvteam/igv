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
import org.apache.commons.lang.StringUtils;

import java.io.File;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Write SAM/BAM Alignments to a file or stream
 * <p/>
 * User: jacob
 * Date: 2012/05/04
 */
public class SAMWriter {

    private static final String SAM_FIELD_SEPARATOR = "\t";

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
            writer = new BAMFileWriter(stream, null);
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

    private static int getFlags(Alignment alignment) {
        int result = alignment.isPaired() ? 0x1 : 0;
        ReadMate mate = alignment.getMate();
        if (mate != null) {
            result += !mate.isMapped() ? 0x8 : 0;
            result += mate.isNegativeStrand() ? 0x20 : 0;
        }
        result += alignment.isProperPair() ? 0x2 : 0;
        result += !alignment.isMapped() ? 0x4 : 0;
        result += alignment.isNegativeStrand() ? 0x10 : 0;
        result += alignment.isFirstOfPair() ? 0x40 : 0;
        result += alignment.isSecondOfPair() ? 0x80 : 0;
        //TODO Not really clear on the meaning of this flag : it seems like we
        //can do without it though
        //result += false ? 0x100 : 0;
        result += alignment.isVendorFailedRead() ? 0x200 : 0;
        result += alignment.isDuplicate() ? 0x400 : 0;
        return result;
    }

    /**
     * Create SAM string from alignment. Work in progress.
     * Currently ignores the quality string and any optional attributes,
     * but should otherwise be correct.
     */
    public static String getSAMString(Alignment alignment) {

        String refName = alignment.getChr();
        List<String> tokens = new ArrayList<String>(11);

        tokens.add(alignment.getReadName());
        tokens.add(Integer.toString(getFlags(alignment)));
        tokens.add(refName);
        tokens.add(Integer.toString(alignment.getAlignmentStart()));
        tokens.add(Integer.toString(alignment.getMappingQuality()));
        tokens.add(alignment.getCigarString());

        ReadMate mate = alignment.getMate();
        String mateRefName = mate != null ? mate.getChr() : null;
        if (refName.equals(mateRefName) &&
                !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(mateRefName)) {
            tokens.add("=");
        } else {
            tokens.add(mateRefName);
        }

        int mateStart = mate != null ? mate.getStart() : 0;
        tokens.add(Integer.toString(mateStart));
        tokens.add(Integer.toString(alignment.getInferredInsertSize()));
        tokens.add(alignment.getReadSequence());

        //TODO Implement quality
        tokens.add("*");
        //tokens.add(SAMUtils.phredToFastq(alignment.getQualityArray()));

        //We add a newline to be consistent with samtools
        String out = StringUtils.join(tokens, SAM_FIELD_SEPARATOR) + "\n";
        return out;

        //TODO Most of our alignment implementations don't have these attributes
//            SAMBinaryTagAndValue attribute = alignment.getBinaryAttributes();
//            while (attribute != null) {
//                out.write(FIELD_SEPARATOR);
//                final String encodedTag;
//                if (attribute.isUnsignedArray()) {
//                    encodedTag = tagCodec.encodeUnsignedArray(tagUtil.makeStringTag(attribute.tag), attribute.value);
//                } else {
//                    encodedTag = tagCodec.encode(tagUtil.makeStringTag(attribute.tag), attribute.value);
//                }
//                out.write(encodedTag);
//                attribute = attribute.getNext();
//            }

    }

    /**
     * Takes an iterable of Alignments, and returns an iterable
     * consisting only of the SamAlignments contained therein.
     */
    public static class SamAlignmentIterable implements Iterable<SamAlignment>, Iterator<SamAlignment> {

        private Iterator<Alignment> alignments;
        private SamAlignment nextAlignment;

        public SamAlignmentIterable(Iterable<Alignment> alignments) {
            this.alignments = alignments.iterator();
            advance();
        }

        private void advance() {
            Alignment next;
            nextAlignment = null;
            while (alignments.hasNext() && nextAlignment == null) {
                next = alignments.next();
                if (next instanceof SamAlignment) {
                    nextAlignment = (SamAlignment) next;
                }
            }
        }

        @Override
        public boolean hasNext() {
            return nextAlignment != null;
        }

        @Override
        public SamAlignment next() {
            if(!hasNext()) throw new NoSuchElementException("No more SamAlignments");
            SamAlignment next = nextAlignment;
            advance();
            return next;
        }

        @Override
        public void remove() {
            //pass
        }

        @Override
        public Iterator iterator() {
            return this;
        }
    }


}
