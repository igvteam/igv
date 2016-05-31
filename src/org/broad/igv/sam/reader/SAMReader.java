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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.EmptyAlignmentIterator;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * A wrapper for SamTextReader that supports query by interval.
 *
 * @author jrobinso
 */
public class SAMReader implements AlignmentReader<PicardAlignment> {

    static Logger log = Logger.getLogger(SAMReader.class);
    String samFile;
    FeatureIndex featureIndex;
    SAMFileHeader header;
    List<String> sequenceNames;

    public SAMReader(String samFile) throws IOException {
        this(samFile, true);
    }

    public SAMReader(String samFile, boolean requireIndex) throws IOException {
        this.samFile = samFile;
        loadHeader();

        if (requireIndex) {
            featureIndex = SamUtils.getIndexFor(samFile);
            if (featureIndex == null) {
                throw new IndexNotFoundException(samFile);
            }
        }
    }

    public SAMFileHeader getFileHeader() {
        return header;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getFileHeader());
    }

    private void loadHeader() {

        InputStream is = null;
        SAMFileReader reader = null;
        try {
            is = ParsingUtils.openInputStreamGZ(new ResourceLocator(samFile));
            BufferedInputStream bis = new BufferedInputStream(is);
            SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
            reader = new SAMFileReader(bis);
            header = reader.getFileHeader();
        } catch (IOException e) {
            log.error("Error loading header", e);
        } finally {
            try {
                if (is != null) {
                    is.close();
                }
                if (reader != null) {
                    reader.close();
                }

            } catch (Exception e) {

            }
        }
    }

    public CloseableIterator<PicardAlignment> query(final String sequence, final int start, final int end, final boolean contained) {

        if (featureIndex == null) {
            featureIndex = SamUtils.getIndexFor(samFile);
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
        int startTileNumber = Math.max(0, (start - startAdjustment)) / featureIndex.getTileWidth();

        FeatureIndex.TileDef seekPos = featureIndex.getTileDef(sequence, startTileNumber);

        if (seekPos != null) {
            // Skip to the start of the query interval and open a sam file reader
            SAMFileReader reader = getSAMFileReader(samFile, seekPos.getStartPosition());
            CloseableIterator<SAMRecord> iter = reader.iterator();
            return new SAMQueryIterator(sequence, start, end, contained, iter);
        }
        return EmptyAlignmentIterator.getInstance();
    }

    public boolean hasIndex() {
        if (featureIndex == null) {
            getIndex();
        }
        return featureIndex != null;
    }

    public void close() throws IOException {
        // Nothing to close
    }


    private FeatureIndex getIndex() {
        if (featureIndex == null) {
            featureIndex = SamUtils.getIndexFor(samFile);
        }
        return featureIndex;
    }

    public List<String> getSequenceNames() {
        if (sequenceNames == null) {
            FeatureIndex idx = getIndex();
            if (idx == null) {
                return null;
            } else {
                sequenceNames = new ArrayList<String>(idx.getIndexedChromosomes());
            }
        }
        return sequenceNames;

    }

    public CloseableIterator<PicardAlignment> iterator() {
        SAMFileReader reader = getSAMFileReader(samFile, -1);
        CloseableIterator<SAMRecord> iter = reader.iterator();
        return new SAMQueryIterator(iter);
    }

    private SAMFileReader getSAMFileReader(String samFile, long startPosition) {
        try {
            SeekableStream stream = IGVSeekableStreamFactory.getInstance().getStreamFor(samFile);
            if (startPosition >= 0) {
                stream.seek(startPosition);
            }
            SAMFileReader reader = new SAMFileReader(stream);
            reader.setValidationStringency(ValidationStringency.SILENT);

            //Need to keep the file source, if loading lazily
            //TODO Can't reload from SAM files. See SAMTextReader.getIterator
            //reader.enableFileSource(PicardAlignment.DEFAULT_LAZY_LOAD);

            return reader;
        } catch (IOException ex) {
            log.error("Error opening sam file", ex);
            throw new RuntimeException("Error opening: " + samFile, ex);
        }
    }

}
