/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import net.sf.samtools.SAMFileHeader;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.tribble.readers.AsciiLineReader;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.EmptyAlignmentIterator;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Set;

/**
 * A wrapper for SamTextReader that supports query by interval.
 *
 * @author jrobinso
 */
public class GeraldQueryReader implements AlignmentQueryReader {

    private static Logger log = Logger.getLogger(GeraldQueryReader.class);
    static int MAX_READ_LENGTH = 100;
    static int maxTileCount = 20;
    String alignmentFile;
    FeatureIndex featureIndex;
    FileInputStream is;
    AlignmentParser parser;


    public GeraldQueryReader(String alignmentFile) {
        this(alignmentFile, true);
    }

    public GeraldQueryReader(String alignmentFile, boolean requireIndex) {
        this.alignmentFile = alignmentFile;
        parser = getParserFor(alignmentFile);
        try {
            is = new FileInputStream(alignmentFile);
        } catch (FileNotFoundException fileNotFoundException) {
            fileNotFoundException.printStackTrace();
        }
        if (requireIndex) {
            featureIndex = SamUtils.getIndexFor(alignmentFile);
            if (featureIndex == null) {
                throw new DataLoadException("Could not locate index file.", alignmentFile);
            }
        }
    }

    // 

    private static AlignmentParser getParserFor(String fn) {
        if (fn.endsWith(".aligned") || (fn.endsWith(".aligned.txt"))) {
            return new DotAlignedParser();
        } else if (fn.endsWith(".bedz") || fn.endsWith(".bed")) {
            return new DotAlignedParser(true);
        } else if (fn.endsWith(".psl") || fn.endsWith(".psxl")) {
            return new PSLAlignmentParser();
        } else {
            return new GeraldParser();
        }
    }

    public SAMFileHeader getHeader() {
        return null;
    }


    private FeatureIndex getIndex() {
        if (featureIndex == null) {
            featureIndex = SamUtils.getIndexFor(alignmentFile);
        }
        return featureIndex;
    }

    public Set<String> getSequenceNames() {
        FeatureIndex idx = getIndex();
        if (idx == null) {
            return null;
        } else {
            return idx.getIndexedChromosomes();
        }
    }

    public CloseableIterator<Alignment> query(final String sequence, final int start, final int end, final boolean contained) {

        if (featureIndex == null) {
            featureIndex = SamUtils.getIndexFor(alignmentFile);
        }

        if (featureIndex == null) {
            throw new java.lang.UnsupportedOperationException("SAM files must be indexed to support query methods");
        }
        if (!featureIndex.containsChromosome(sequence)) {
            return EmptyAlignmentIterator.getInstance();
        }

        // If contained == false (include overlaps) we need to adjust the start to
        // ensure we get features that extend into this segment.
        int startAdjustment = contained ? 0 : featureIndex.getLongestFeature(sequence);
        return new GeraldQueryIterator(sequence, start - startAdjustment, end, contained);

    }

    // TODO -- implementation

    public boolean hasIndex() {
        return true;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void close() throws IOException {
        if (is != null) {
            is.close();
        }
    }

    public CloseableIterator<Alignment> iterator() {
        return new GeraldQueryIterator();
    }

    /**
     *
     */
    class GeraldQueryIterator implements CloseableIterator<Alignment> {

        String chr;
        int start;
        int end;
        boolean contained;
        Alignment currentRecord;
        AsciiLineReader reader;

        public GeraldQueryIterator() {
            this.chr = null;
            reader = new AsciiLineReader(is);
            readNextRecord();
        }

        public GeraldQueryIterator(String sequence, int start, int end, boolean contained) {
            this.chr = sequence;
            this.start = start;
            this.end = end;
            this.contained = contained;
            seekToStart();
            reader = new AsciiLineReader(is);
            advanceToFirstRecord();
        }

        private Alignment readNextRecord() {
            try {
                currentRecord = parser.readNextRecord(reader);
            } catch (IOException e) {
                log.error("Error reading Gerald record", e);
                return null;
            }
            return currentRecord;
        }

        private void advanceToFirstRecord() {
            while ((readNextRecord()) != null) {
                if (!currentRecord.getChr().equals(chr)) {
                    break;
                } else if ((contained && currentRecord.getStart() >= start) ||
                        (!contained && currentRecord.getEnd() >= start)) {
                    break;
                }
            }
        }

        public void close() {
            try {
                is.close();
            } catch (IOException ex) {
                log.error("Error closing alignment file", ex);
            }
        }

        public boolean hasNext() {

            // TODO chr == null implies an iterator, not query.  Fix this
            if (chr == null) {
                return currentRecord != null;
            }

            if (currentRecord == null ||
                    !chr.equals(currentRecord.getChr())) {
                return false;
            } else {
                return contained ? currentRecord.getEnd() <= end
                        : currentRecord.getStart() <= end;
            }
        }

        public Alignment next() {
            Alignment ret = currentRecord;
            readNextRecord();
            return ret;

        }

        public void remove() {
            throw new UnsupportedOperationException("Not supported yet.");
        }

        private void seekToStart() {

            if (featureIndex == null) {
                throw new java.lang.UnsupportedOperationException("SAM files must be indexed to support query methods");
            }

            // If contained == false (include overlaps) we need to adjust the start to
            // ensure we get features that extend into this segment.
            int startAdjustment = contained ? 0 : MAX_READ_LENGTH; //featureIndex.getLongestFeature(chr);
            int startTileNumber = Math.max(0, (start - startAdjustment)) / featureIndex.getTileWidth();

            FeatureIndex.TileDef seekPos = featureIndex.getTileDef(chr, startTileNumber);
            long startPosition = seekPos == null ? 0 : seekPos.getStartPosition();

            try {
                // Skip to the start of the chromosome (approximate)
                is = new FileInputStream(alignmentFile);
                is.getChannel().position(startPosition);

            } catch (Exception ex) {
                throw new RuntimeException(ex);
            }
        }
    }

}
