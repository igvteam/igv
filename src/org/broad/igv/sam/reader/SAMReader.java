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
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.EmptyAlignmentIterator;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

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
public class SAMReader implements AlignmentReader {

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

    public SAMFileHeader getHeader() {
        return header;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getHeader());
    }

    private void loadHeader() {

        InputStream is = null;
        SAMFileReader reader = null;
        try {
            is = ParsingUtils.openInputStreamGZ(new ResourceLocator(samFile));
            BufferedInputStream bis = new BufferedInputStream(is);
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

    public CloseableIterator<Alignment> query(final String sequence, final int start, final int end, final boolean contained) {

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

    public CloseableIterator<Alignment> iterator() {
        SAMFileReader reader = getSAMFileReader(samFile, -1);
        CloseableIterator<SAMRecord> iter = reader.iterator();
        return new SAMQueryIterator(iter);
    }

    private SAMFileReader getSAMFileReader(String samFile, long startPosition) {
        try {
            SeekableStream stream = SeekableStreamFactory.getStreamFor(samFile);
            if (startPosition >= 0) {
                stream.seek(startPosition);
            }
            SAMFileReader reader = new SAMFileReader(stream);
            reader.setValidationStringency(ValidationStringency.SILENT);

            //Need to keep the file source, if loading lazily
            //TODO Can't reload from SAM files. See SAMTextReader.getIterator
            //reader.enableFileSource(SamAlignment.DEFAULT_LAZY_LOAD);

            return reader;
        } catch (IOException ex) {
            log.error("Error opening sam file", ex);
            throw new RuntimeException("Error opening: " + samFile, ex);
        }
    }

}
