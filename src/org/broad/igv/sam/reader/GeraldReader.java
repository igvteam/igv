/*
* Copyright (c) 2007-2012 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
*
* This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
* Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
*
* THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
* WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
* WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
* PURPOSE, NON-INFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
* OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
* TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
* OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
* ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
* THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
* SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*/

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.EmptyAlignmentIterator;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * A wrapper for SamTextReader that supports query by interval.
 *
 * @author jrobinso
 */
public class GeraldReader implements AlignmentReader {

    private static Logger log = Logger.getLogger(GeraldReader.class);
    static int MAX_READ_LENGTH = 100;
    static int maxTileCount = 20;
    String alignmentFile;
    FeatureIndex featureIndex;
    FileInputStream is;
    AlignmentParser parser;

    public GeraldReader(String alignmentFile, boolean requireIndex) {
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

    private FeatureIndex getIndex() {
        if (featureIndex == null) {
            featureIndex = SamUtils.getIndexFor(alignmentFile);
        }
        return featureIndex;
    }

    public List<String> getSequenceNames() {
        FeatureIndex idx = getIndex();
        if (idx == null) {
            return null;
        } else {
            return new ArrayList<String>(idx.getIndexedChromosomes());
        }
    }

    @Override
    public Set<String> getPlatforms() {
        return null;
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

        return new GeraldQueryIterator(sequence, start, end, contained);

    }

    public boolean hasIndex() {
        return getIndex() != null;
    }

    public void close() throws IOException {
        if (is != null) {
            is.close();
        }
    }

    public CloseableIterator<Alignment> iterator() {
        return new GeraldIterator();
    }

    /**
     * Class which doesn't support any querying, just iterates
     * over all records.
     */
    class GeraldIterator implements CloseableIterator<Alignment> {

        String chr;
        int start;
        int end;
        boolean contained;
        Alignment nextRecord;
        AsciiLineReader reader;

        public GeraldIterator() {
            this.chr = null;
            reader = new AsciiLineReader(is);
            readNextRecord();
        }

        protected Alignment parseNextRecord() {
            Alignment next;
            try {
                next = parser.readNextRecord(reader);
            } catch (IOException e) {
                log.error("Error reading Gerald record", e);
                return null;
            }
            return next;
        }

        protected Alignment readNextRecord() {
            nextRecord = parseNextRecord();
            return nextRecord;
        }

        public void close() {
            try {
                is.close();
            } catch (IOException ex) {
                log.error("Error closing alignment file", ex);
            }
        }

        @Override
        public boolean hasNext() {
            return nextRecord != null;
        }

        @Override
        public Alignment next() {
            Alignment ret = nextRecord;
            readNextRecord();
            return ret;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Not supported yet.");
        }
    }

    /**
     * Iterator for querying.
     * author Jacob Silterra
     */
    class GeraldQueryIterator extends GeraldIterator {

        public GeraldQueryIterator(String sequence, int start, int end, boolean contained) {
            this.chr = sequence;
            this.start = start;
            this.end = end;
            this.contained = contained;
            seekToStart();
            reader = new AsciiLineReader(is);
            advanceToFirstRecord();
        }

        private boolean withinBounds(Alignment alignment) {
            boolean within = alignment.getEnd() <= end && alignment.getStart() >= start;
            if (contained) {
                return within;
            }
            within |= alignment.getStart() <= start && alignment.getEnd() > start;
            within |= alignment.getStart() >= start && alignment.getStart() < end;
            return within;
        }

        private void advanceToFirstRecord() {
            nextRecord = parseNextRecord();
            while (nextRecord != null) {
                if (!nextRecord.getChr().equals(chr)) {
                    break;
                } else if (withinBounds(nextRecord)) {
                    break;
                }
                readNextRecord();
            }
        }

        @Override
        protected Alignment readNextRecord() {
            advance();
            //We use exclusive end
            while ((nextRecord != null) && !withinBounds(nextRecord)) {
                advance();
            }
            return nextRecord;
        }

        private void advance() {
            if (hasNext()) {
                nextRecord = parseNextRecord();
                if (nextRecord == null) {
                    return;
                }

                if (nextRecord.getStart() >= end) {
                    nextRecord = null;
                }
            } else {
                nextRecord = null;
            }
        }

        @Override
        public boolean hasNext() {
            if (nextRecord == null ||
                    !chr.equals(nextRecord.getChr())) {
                return false;
            } else {
                return contained ? nextRecord.getEnd() <= end
                        : nextRecord.getStart() < end;
            }
        }

        private void seekToStart() {

            if (featureIndex == null) {
                throw new java.lang.UnsupportedOperationException("SAM files must be indexed to support query methods");
            }

            // If contained == false (include overlaps) we need to adjust the start to
            // ensure we get features that extend into this segment.
            int startAdjustment = contained ? 0 : featureIndex.getLongestFeature(chr);
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
