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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;
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
