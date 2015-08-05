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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.EmptyAlignmentIterator;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.ui.util.MessageUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 */
public class BAMFileReader implements AlignmentReader<PicardAlignment> {

    private static Logger log = Logger.getLogger(BAMFileReader.class);
    SAMFileReader reader;
    SAMFileHeader header;

    public BAMFileReader(File bamFile) {
        try {
            File indexFile = findIndexFile(bamFile);
            reader = new SAMFileReader(bamFile, indexFile);
            reader.setValidationStringency(ValidationStringency.SILENT);
            loadHeader();
        } catch (Exception e) {
            MessageUtils.showMessage("Error loading SAM header: " + e.getMessage());
            e.printStackTrace();
        }
    }

    public SAMFileHeader getFileHeader() {
        if (header == null) {
            loadHeader();
        }
        return header;
    }

    public List<String> getSequenceNames() {
        SAMFileHeader header = getFileHeader();
        if (header == null) {
            return null;
        }
        List<String> seqNames = new ArrayList();
        List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
        if (records.size() > 0) {
            for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                String chr = rec.getSequenceName();
                seqNames.add(chr);
            }
        }
        return seqNames;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getFileHeader());
    }

    private void loadHeader() {
        header = reader.getFileHeader();
    }

    public void close() throws IOException {
        reader.close();
    }

    public boolean hasIndex() {
        return reader.hasIndex();
    }

    public CloseableIterator<PicardAlignment> query(String sequence, int start, int end, boolean contained) {
        SAMRecordIterator query = null;
        try {
            query = reader.query(sequence, start + 1, end, contained);
            return new WrappedIterator(query);
        } catch (ArrayIndexOutOfBoundsException e) {
            log.error("Error querying BAM file ", e);
            MessageUtils.showMessage("Error reading bam file.  This usually indicates a problem with the index (bai) file." +
                    "<br>" + e.toString() + " (" + e.getMessage() + ")");
            return EmptyAlignmentIterator.getInstance();
        }

    }

    public CloseableIterator<PicardAlignment> iterator() {
        return new WrappedIterator(reader.iterator());
    }

    /**
     * Look for BAM index file according to standard naming convention.  Slightly modified version of Picard
     * function of the same name.
     *
     * @param dataFile BAM file name.
     * @return Index file name, or null if not found.
     */
    private static File findIndexFile(final File dataFile) {

        final String bamPath = dataFile.getAbsolutePath();

        // foo.bam.bai
        String bai = bamPath + ".bai";
        File indexFile1 = new File(bai);
        if (indexFile1.length() > 0) {
            return indexFile1;
        }

        // alternate (Picard) convention,  foo.bai
        final String bamExtension = ".bam";
        File indexFile2 = null;
        if (bamPath.toLowerCase().endsWith(bamExtension)) {
            bai = bamPath.substring(0, bamPath.length() - bamExtension.length()) + ".bai";
            indexFile2 = new File(bai);
            if (indexFile2.length() > 0) {
                return indexFile2;
            }
        }

        log.info("Index file: " + indexFile1.getAbsolutePath() + " not found");
        if (indexFile2 != null) log.info("Index file: " + indexFile2.getAbsolutePath() + " not Found");

        return null;
    }


}
