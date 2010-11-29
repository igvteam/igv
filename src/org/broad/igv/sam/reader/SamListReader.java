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

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.*;

/**
 * Poorly named class.  Created to wrap a list of sam files (1 per chromosome) as a single logical file
 */

public class SamListReader implements AlignmentQueryReader {

    private static Logger log = Logger.getLogger(SamListReader.class);

    SAMFileHeader header;

    Map<String, ResourceLocator> locators;
    Map<String, AlignmentQueryReader> readers;

    public SamListReader(Map<String, ResourceLocator> locators) {
        this.locators = locators;
        readers = new HashMap(locators.size());
    }

    public void close() throws IOException {
        for (AlignmentQueryReader reader : readers.values()) {
            reader.close();
        }
        readers.clear();
    }

    public SAMFileHeader getHeader() throws IOException {
        if (locators == null || locators.size() == 0) {
            return null;
        }
        if (header == null) {
            // Any header will do, see if we have one
            if (readers.size() == 0) {
                Map.Entry<String, ResourceLocator> tmp = locators.entrySet().iterator().next();
                AlignmentQueryReader firstReader = SamQueryReaderFactory.getReader(tmp.getValue());
                readers.put(tmp.getKey(), firstReader);
            }
            header = (readers.values().iterator().next()).getHeader();
        }
        return header;
    }

    public Set<String> getSequenceNames() {
        try {
            SAMFileHeader header = getHeader();
            if (header == null) {
                return null;
            }
            Set<String> seqNames = new HashSet();
            List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
            if (records.size() > 0) {
                for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                    String chr = rec.getSequenceName();
                    seqNames.add(chr);
                }
            }
            return seqNames;
        } catch (IOException e) {
            log.error("Error fetching sequence names", e);
            return null;
        }
    }

    public CloseableIterator<Alignment> iterator() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public CloseableIterator<Alignment> query(String chr, int start, int end, boolean contained) throws IOException {
        AlignmentQueryReader reader = readers.get(chr);
        if (reader == null) {
            ResourceLocator locator = locators.get(chr);
            if (locator == null) {
                // TODO -- return empty iterator?
                return null;
            }
            reader = SamQueryReaderFactory.getReader(locator);
        }
        return reader.query(chr, start, end, contained);
    }

    public boolean hasIndex() {
        // TODO -- test one of the readers
        return true;
    }
}
