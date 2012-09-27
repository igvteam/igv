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

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.util.MessageUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 */
public class BAMFileReader implements AlignmentReader {

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

    private SAMFileHeader getHeader() {
        if (header == null) {
            loadHeader();
        }
        return header;
    }

    public List<String> getSequenceNames() {
        SAMFileHeader header = getHeader();
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
        return AlignmentReaderFactory.getPlatforms(getHeader());
    }

    private void loadHeader() {
        header = reader.getFileHeader();
    }

    public void close() throws IOException {
        reader.close();
    }

    public boolean hasIndex() {
        return reader.getIndex() != null;
    }

    public CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) {
        SAMRecordIterator query = null;
        try {
            query = reader.query(sequence, start + 1, end, contained);
            return new WrappedIterator(query);
        } catch (ArrayIndexOutOfBoundsException e) {
            log.error("Error querying BAM file ", e);
            MessageUtils.showMessage("Error reading bam file.  This usually indicates a problem with the index (bai) file." +
                    "<br>" + e.toString() + " (" + e.getMessage() + ")");
            return EmptyIterator.instance;
        }

    }

    public CloseableIterator<Alignment> iterator() {
        return new WrappedIterator(reader.iterator());
    }


    /**
     * Class
     */
    public static class EmptyIterator implements CloseableIterator<Alignment> {

        public static final EmptyIterator instance = new EmptyIterator();

        public void close() {
            // Ignore
        }

        public boolean hasNext() {
            return false;
        }

        public Alignment next() {
            return null;
        }

        public void remove() {

        }
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
