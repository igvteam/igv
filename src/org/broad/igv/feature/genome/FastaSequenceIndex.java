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
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Set;

/**
 * Representation of a fasta (.fai) index.  This is a modified version of a similar class in Picard.
 * <p/>
 * author: jrobinso
 * Date: 8/7/11
 */
public class FastaSequenceIndex {

    static Logger log = Logger.getLogger(FastaSequenceIndex.class);

    /**
     * Store the entries.  Use a LinkedHashMap for consistent iteration in insertion order.
     */
    private final LinkedHashMap<String, FastaSequenceIndexEntry> sequenceEntries = new LinkedHashMap<String, FastaSequenceIndexEntry>();

    public FastaSequenceIndex(String indexPath) throws IOException {
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
     * Creates an index for the provided fasta file
     * inputPath can be a URL, outputPath must point to a file.
     *
     * @param inputPath
     * @param outputPath
     * @return
     * @throws DataLoadException If the fasta file cannot be indexed, for instance
     *                           because the lines are of an uneven length
     */
    public static void createIndexFile(String inputPath, String outputPath) throws DataLoadException, IOException {

        AsciiLineReader reader = null;
        BufferedWriter writer = null;

        reader = new AsciiLineReader(ParsingUtils.openInputStream(inputPath));
        writer = new BufferedWriter(new FileWriter(outputPath));
        String line = null;
        String curContig = null;
        int basesPerLine = -1, bytesPerLine = -1;
        long location = 0, size = 0, lastPosition = 0;

        int basesThisLine, bytesThisLine;
        int numInconsistentLines = -1;
        boolean haveTasks = true;


        //We loop through, generating a new FastaSequenceIndexEntry
        //every time we see a new header line, or when the file ends.
        while (haveTasks) {
            line = reader.readLine();
            //Treat empty line as end of file
            //This can come up for trailing newline
            if (line == null || line.trim().length() == 0) {
                line = null;
            }
            if (line == null || line.startsWith(">")) {
                //The last line can have a different number of bases/bytes
                if (numInconsistentLines >= 2) {
                    throw new DataLoadException("Fasta file has uneven line lengths in contig " + curContig, inputPath);
                }

                //Done with old contig
                if (curContig != null) {
                    writeLine(writer, curContig, size, location, basesPerLine, bytesPerLine);
                }

                if (line == null) {
                    haveTasks = false;
                    break;
                }

                //Header line
                curContig = line.split("\\s")[0];
                curContig = curContig.substring(1);
                //Should be starting position of next line
                location = reader.getPosition();
                size = 0;
                basesPerLine = -1;
                bytesPerLine = -1;
                numInconsistentLines = -1;
            } else {
                basesThisLine = line.length();
                bytesThisLine = (int) (reader.getPosition() - lastPosition);

                //Calculate stats per line if first line, otherwise
                //check for consistency
                if (numInconsistentLines < 0) {
                    basesPerLine = basesThisLine;
                    bytesPerLine = bytesThisLine;
                    numInconsistentLines = 0;
                } else {
                    if (basesPerLine != basesThisLine || bytesPerLine != bytesThisLine) {
                        numInconsistentLines++;
                    }
                }

                size += basesThisLine;
            }
            lastPosition = reader.getPosition();
        }

        reader.close();
        writer.close();
    }

    private static void writeLine(Writer writer, String contig, long size, long location, int basesPerLine, int bytesPerLine) throws IOException {
        String delim = "\t";
        String line = contig + delim + size + delim + location + delim + basesPerLine + delim + bytesPerLine;
        writer.write(line);
        //We infer the newline character based on bytesPerLine - basesPerLine
        //Fasta file may not have been created on this platform, want to keep the index and fasta file consistent
        String newline = "\n";
        if (bytesPerLine - basesPerLine == 2) {
            newline = "\r\n";
        }
        writer.write(newline);
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
        String[] tokens = new String[5];
        try {
            reader = ParsingUtils.openBufferedReader(indexFile);

            String nextLine;
            while ((nextLine = reader.readLine()) != null) {

                // Tokenize and validate the index line.
                int nTokens = ParsingUtils.splitWhitespace(nextLine, tokens);
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