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

package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Set;

/**
 * Representation of a fasta (.fai) index.  This is a modified version of a similar class in Picard, but extended
 * to handle loading by URLs.
 *
 * author: jrobinso
 * Date: 8/7/11
 */
public class FastaIndex {

    static Logger log = Logger.getLogger(FastaIndex.class);

    /**
     * Store the entries.  Use a LinkedHashMap for consistent iteration in insertion order.
     */
    private final LinkedHashMap<String, FastaSequenceIndexEntry> sequenceEntries = new LinkedHashMap<String, FastaSequenceIndexEntry>();

    public FastaIndex(String indexPath) throws IOException {
        parseIndexFile(indexPath);
    }

    public Set<String> getSequenceNames() {
        return sequenceEntries.keySet();
    }

    public FastaSequenceIndexEntry getIndexEntry(String name) {
        return sequenceEntries.get(name);
    }


    public int getSequenceSize(String name) {
        FastaSequenceIndexEntry entry = sequenceEntries.get(name);
        return entry == null ? -1 : (int) entry.getSize();
    }

    /**
     * Parse the contents of an index file
     * <p/>
     * Example index file
     * <p/>
     * chr01p	6220112	8	50	51
     * chr02q	8059593	6344531	50	51
     * chr03q	5803340	14565324	50	51
     *
     * @param indexFile File to parse.
     * @throws java.io.FileNotFoundException Thrown if file could not be opened.
     */
    private void parseIndexFile(String indexFile) throws IOException {

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(indexFile);

            String nextLine;
            while ((nextLine = reader.readLine()) != null) {

                // Tokenize and validate the index line.
                String[] tokens =  Globals.singleTabMultiSpacePattern.split(nextLine);
                int nTokens =  tokens.length;
                if (nTokens != 5) {
                    throw new RuntimeException("Error.  Unexpected number of tokens parsing: " + indexFile);
                }
                // Parse the index line.
                String contig = tokens[0];
                contig = GenomeImporter.SEQUENCE_NAME_SPLITTER.split(contig, 2)[0];
                long size = Long.parseLong(tokens[1]);
                long location = Long.parseLong(tokens[2]);
                int basesPerLine = Integer.parseInt(tokens[3]);
                int bytesPerLine = Integer.parseInt(tokens[4]);

                // Build sequence structure
                add(new FastaSequenceIndexEntry(contig, location, size, basesPerLine, bytesPerLine));
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    private void add(FastaSequenceIndexEntry indexEntry) {
        final FastaSequenceIndexEntry ret = sequenceEntries.put(indexEntry.getContig(), indexEntry);
        if (ret != null) {
            throw new RuntimeException("Contig '" + indexEntry.getContig() + "' already exists in fasta index.");
        }
    }


    /**
     * Hold an individual entry in a fasta sequence index file.
     */
    static class FastaSequenceIndexEntry {
        private String contig;
        private long position;
        private long size;
        private int basesPerLine;
        private int bytesPerLine;


        /**
         * Create a new entry with the given parameters.
         *
         * @param contig       Contig this entry represents.
         * @param position     Location (byte coordinate) in the fasta file.
         * @param size         The number of bases in the contig.
         * @param basesPerLine How many bases are on each line.
         * @param bytesPerLine How many bytes are on each line (includes newline characters).
         */
        public FastaSequenceIndexEntry(String contig,
                                       long position,
                                       long size,
                                       int basesPerLine,
                                       int bytesPerLine) {
            this.contig = contig;
            this.position = position;
            this.size = size;
            this.basesPerLine = basesPerLine;
            this.bytesPerLine = bytesPerLine;
        }

        /**
         * Gets the contig associated with this entry.
         *
         * @return String representation of the contig.
         */
        public String getContig() {
            return contig;
        }

        /**
         * Gets the position of this contig within the fasta.
         *
         * @return seek position within the fasta.
         */
        public long getPosition() {
            return position;
        }

        /**
         * Gets the size, in bytes, of the data in the contig.
         *
         * @return size of the contig bases in bytes.
         */
        public long getSize() {
            return size;
        }

        /**
         * Gets the number of bases in a given line.
         *
         * @return Number of bases in the fasta line.
         */
        public int getBasesPerLine() {
            return basesPerLine;
        }

        /**
         * How many bytes (bases + whitespace) are consumed by the
         * given line?
         *
         * @return Number of bytes in a line.
         */
        public int getBytesPerLine() {
            return bytesPerLine;
        }


        /**
         * For debugging.  Emit the contents of each contig line.
         *
         * @return A string representation of the contig line.
         */
        public String toString() {
            return String.format("contig %s; position %d; size %d; basesPerLine %d; bytesPerLine %d", contig,
                    position,
                    size,
                    basesPerLine,
                    bytesPerLine);
        }

    }
}