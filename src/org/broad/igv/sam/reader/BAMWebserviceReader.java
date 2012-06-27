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
import net.sf.samtools.util.BufferedLineReader;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.LineReader;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SamAlignment;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/**
 * Reader for querying a "bam" webservice.
 *
 * @author jrobinso
 */
public class BAMWebserviceReader implements AlignmentReader {

    Logger log = Logger.getLogger(BAMWebserviceReader.class);
    String serverURL;
    String file;
    SAMFileHeader header;
    List<String> sequenceNames;

    public BAMWebserviceReader(ResourceLocator locator) {
        this.serverURL = locator.getServerURL();
        this.file = locator.getPath();
        loadHeader();
    }

    public void close() throws IOException {
        // Nothing to do
    }

    public boolean hasIndex() {
        // There is no server API to to check this, so we'll be optimisitic.
        return true;
    }

    public CloseableIterator<Alignment> query(String chr, int start, int end, boolean contained) {
        try {
            URL url = new URL(serverURL + "?method=samQuery&samFile=" + file + "&chr=" +
                    chr + "&start=" + start + "&end=" + end + "&contained=" + contained);
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);
            return new RemoteQueryIterator(new GZIPInputStream(is, 8192));

        } catch (IOException ex) {
            log.error("Error opening file", ex);
            throw new RuntimeException(ex);
        }
    }

    public SAMFileHeader getHeader() {
        if (header == null) {
            loadHeader();
        }
        return header;
    }

    public List<String> getSequenceNames() {
        if (sequenceNames == null) {
            SAMFileHeader header = getHeader();
            if (header == null) {
                return null;
            }
            sequenceNames = new ArrayList();
            List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
            if (records.size() > 0) {
                boolean ensembleChrConventions = true;
                for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                    String chr = rec.getSequenceName();
                    sequenceNames.add(chr);
                }

            }
        }
        return sequenceNames;
    }


    /**
     * @return true if any readgroups have the platform tag set to "IONTORRENT"
     */
    public Set<String> getPlatforms() {
        Set<String> platforms = new HashSet<String>();
        SAMFileHeader header = getHeader();
        if (header != null) {
            List<SAMReadGroupRecord> readGroups = header.getReadGroups();
            if (readGroups != null) {
                platforms = new HashSet<String>();
                for (SAMReadGroupRecord rg : readGroups) {
                    platforms.add(rg.getPlatform());
                }
            }
        }
        return platforms;
    }

    private void loadHeader() {
        InputStream is = null;
        try {
            URL url = new URL(serverURL + "?method=samHeader&samFile=" + file);
            is = HttpUtils.getInstance().openConnectionStream(url);

            LineReader reader = new BufferedLineReader(is);
            SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            header = codec.decode(reader, null);

        } catch (IOException ex) {
            log.error("Error opening file", ex);
            throw new RuntimeException(ex);
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException ex) {
                    log.error("Error closing url stream", ex);
                }
            }
        }
    }

    public CloseableIterator<Alignment> iterator() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    class RemoteQueryIterator implements CloseableIterator<Alignment> {

        InputStream inputStream;
        SAMRecord nextRecord;
        BAMRecordCodec codec;

        public RemoteQueryIterator(InputStream is) {
            this.inputStream = is;
            codec = new BAMRecordCodec(header);
            codec.setInputStream(is);
            advance();
        }

        private void advance() {
            nextRecord = codec.decode();
        }

        public void close() {
            if (inputStream != null) {
                try {
                    inputStream.close();
                    inputStream = null;
                } catch (IOException ex) {
                    log.error("Error closing sam record stream", ex);
                }
            }
        }

        public boolean hasNext() {
            return nextRecord != null;
        }

        public SamAlignment next() {
            SamAlignment ret = new SamAlignment(nextRecord);
            advance();
            return ret;
        }

        public void remove() {
            // ignored
        }

        // Just in case caller forgot to close the iterator

        @Override
        protected void finalize() throws Throwable {
            super.finalize();
            if (inputStream != null) {
                inputStream.close();
            }
        }
    }
}
