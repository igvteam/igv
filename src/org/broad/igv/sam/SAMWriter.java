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

package org.broad.igv.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.feature.Range;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.Utilities;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Write SAM/BAM Alignments to a file or stream
 * <p/>
 * @author jacob
 * @since 2012/05/04
 */
public class SAMWriter {

    private static final String SAM_FIELD_SEPARATOR = "\t";

    private SAMFileHeader header;

    public SAMWriter(SAMFileHeader header) {
        this.header = header;
    }

    public int writeToFile(File outFile, Iterator<PicardAlignment> alignments, boolean createIndex) {
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        factory.setCreateIndex(createIndex);
        SAMFileWriter writer = factory.makeSAMOrBAMWriter(header, true, outFile);
        return writeAlignments(writer, alignments);
    }

    public int writeToStream(OutputStream stream, Iterator<PicardAlignment> alignments, boolean bam) {

        SAMFileWriterImpl writer;
       // if (bam) {
       //     return 0;   // Don't know how to output bams
       //     //writer = new BAMFileWriter(stream, null);
       // } else {
            writer = new SAMTextWriter(stream);

            writer.setHeader(header);
            return writeAlignments(writer, alignments);
       // }
    }

    private int writeAlignments(SAMFileWriter writer, Iterator<PicardAlignment> alignments) {
        int count = 0;
        while (alignments.hasNext()) {
            PicardAlignment al = alignments.next();
            writer.addAlignment(al.getRecord());
            count++;
        }
        writer.close();
        return count;
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
     * Takes an iterator of Alignments, and returns an iterable/iterator
     * consisting only of the SamAlignments contained therein.
     * Can also be used to filter by position
     */
    public static class SamAlignmentIterable implements Iterable<PicardAlignment>, Iterator<PicardAlignment> {

        private Iterator<Alignment> alignments;
        private PicardAlignment nextAlignment;
        private String chr = null;
        private int start = -1;
        private int end = -1;

        public SamAlignmentIterable(Iterator<Alignment> alignments, String chr, int start, int end) {
            this.alignments = alignments;
            this.chr = chr;
            this.start = start;
            this.end = end;
            advance();
        }

        private void advance() {
            Alignment next;
            nextAlignment = null;
            while (alignments.hasNext() && nextAlignment == null) {
                next = alignments.next();
                if (next instanceof PicardAlignment && passLocFilter(next)) {
                    nextAlignment = (PicardAlignment) next;
                }
            }
        }

        @Override
        public boolean hasNext() {
            return nextAlignment != null;
        }

        @Override
        public PicardAlignment next() {
            if(!hasNext()) throw new NoSuchElementException("No more SamAlignments");
            PicardAlignment next = nextAlignment;
            advance();
            return next;
        }

        @Override
        public void remove() {
            //pass
        }

        @Override
        public Iterator<PicardAlignment> iterator() {
            return this;
        }

        private boolean passLocFilter(Alignment al){
            return this.chr != null && this.overlaps(al.getChr(), al.getStart(), al.getEnd());
        }

        /**
         * Determine whether there is any overlap between this interval and the specified interval
         */
        private boolean overlaps(String chr, int start, int end) {
            return Utilities.objectEqual(this.chr, chr) && this.start <= end && this.end >= start;
        }
    }

    /**
     * Use Picard to write alignments which are already stored in memory
     * @param dataManager
     * @param outFile
     * @param sequence
     * @param start
     * @param end
     * @return
     * @throws IOException
     */
    public static int writeAlignmentFilePicard(AlignmentDataManager dataManager, File outFile,
                                               String sequence, int start, int end) throws IOException{

        ResourceLocator inlocator = dataManager.getLocator();
        checkExportableAlignmentFile(inlocator.getTypeString());

        final SAMFileHeader fileHeader = dataManager.getReader().getFileHeader();
        //IGV can only load files sorted in coordinate order, but they aren't always
        //labelled as such.
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        Range range = new Range(sequence, start, end);
        Iterator<Alignment> iter = dataManager.getLoadedInterval(range).getAlignmentIterator();
        Iterator<PicardAlignment> samIter = new SamAlignmentIterable(iter, sequence, start, end);

        SAMWriter writer = new SAMWriter(fileHeader);
        return writer.writeToFile(outFile, samIter, true);
    }

    /**
     * Use Picard to write alignment subset, as read from a file
     * @param inlocator
     * @param outPath
     * @param sequence
     * @param start
     * @param end
     * @return
     */
    public static int writeAlignmentFilePicard(ResourceLocator inlocator, String outPath,
                                                String sequence, int start, int end) throws IOException{

        checkExportableAlignmentFile(inlocator.getTypeString());

        AlignmentReader reader = AlignmentReaderFactory.getReader(inlocator);
        CloseableIterator<PicardAlignment> iter = reader.query(sequence, start, end, false);
        final SAMFileHeader fileHeader = reader.getFileHeader();

        SAMWriter writer = new SAMWriter(fileHeader);
        int count = writer.writeToFile(new File(outPath), iter, true);
        iter.close();

        return count;
    }

    private static void checkExportableAlignmentFile(String typeString){
        String[] validExts = new String[]{".bam", ".sam", ".bam.list", ".sam.list"};
        boolean isValidExt = false;
        for(String validExt: validExts){
            isValidExt |= typeString.endsWith(validExt);
        }
        if(!isValidExt){
            throw new IllegalArgumentException("Input alignment valid not valid for export");
        }
    }

}
