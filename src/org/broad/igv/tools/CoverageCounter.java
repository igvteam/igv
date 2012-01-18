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
package org.broad.igv.tools;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.MergedAlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.ui.filefilters.AlignmentFileFilter;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.stats.Distribution;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;

/**
 * Class to compute coverage on an alignment or feature file.  This class is designed to be instantiated and executed
 * from a single thread.
 */
public class CoverageCounter {

    /**
     * The path to the alignment file being counted.
     */
    private String alignmentFile;

    /**
     * A consumer of data produced by this class,  normally a TDF Preprocessor.
     */
    private DataConsumer consumer;

    /**
     * Window size in base pairs.  Genome is divided into non-overlapping windows of this size.  The counts reported
     * are averages over the window.
     */
    private int windowSize = 1;

    /**
     * Minimum mapping quality.  Alignments with MQ less than this value are filtered.
     */
    private int minMappingQuality = 0;

    /**
     * Strand option.
     * TODO -- not currently used,  implement
     */
    private int strandOption = -1;

    /**
     * Extension factor.  Reads are extended by this amount before counting.   The purpose is to yield an approximate
     * count of fragment "coverage", as opposed to read coverage.  If used, the value should be set to
     * extFactor = averageFragmentSize - averageReadLength
     */
    private int extFactor;

    /**
     * Flag to control treatment of duplicates.  If true duplicates are counted.  The default value is false.
     */
    private boolean includeDuplicates = false;

    private Genome genome;

    /**
     * Optional wig file specifier.  If non-null,  a "wiggle" file is created in addition to the TDF file.
     */
    private File wigFile = null;

    /**
     * Total number of alignments that pass filters and are counted.
     */
    private int totalCount = 0;

    /**
     * The query interval, usually this is null but can be used to restrict the interval of the alignment file that is
     * computed.  The file must be indexed (queryable) if this is not null
     */
    private Locus interval;

    /**
     * Data buffer to pass data to the "consumer" (preprocessor).
     */
    private float[] buffer;


    private boolean computeTDF = true;


    /**
     * @param alignmentFile - path to the file to count
     * @param consumer      - the data consumer, in this case a TDF preprocessor
     * @param windowSize    - window size in bp, counts are performed over this window
     * @param extFactor     - the extension factor, read is artificially extended by this amount
     * @param wigFile       - path to the wig file (optional)
     * @param genome        - the Genome,  used to size chromosomes
     * @param options       - additional coverage options (optional)
     */
    public CoverageCounter(String alignmentFile,
                           DataConsumer consumer,
                           int windowSize,
                           int extFactor,
                           File wigFile,
                           Genome genome,
                           String options) {
        this.alignmentFile = alignmentFile;
        this.consumer = consumer;
        this.windowSize = windowSize;
        this.extFactor = extFactor;
        this.wigFile = wigFile;
        this.genome = genome;
        buffer = strandOption < 0 ? new float[1] : new float[2];

        if (options != null) {
            parseOptions(options);
        }


    }

    /**
     * Parse command-line options.  Perhaps this should be done in the calling program.
     *
     * @param options
     */
    private void parseOptions(String options) {
        String[] opts = options.split(",");
        for (String opt : opts) {
            if (opt.startsWith("d")) {
                includeDuplicates = true;
            } else if (opt.startsWith("m=")) {
                String[] tmp = opt.split("=");
                minMappingQuality = Integer.parseInt(tmp[1]);
                System.out.println("Minimum mapping quality = " + minMappingQuality);
            } else if (opt.startsWith("q")) {
                String[] tmp = opt.split("@");
                interval = new Locus(tmp[1]);
            } else {
                System.out.println("Unknown coverage option: " + opt);
            }
        }

    }


    // TODO -- command-line options to ovveride all of these checks
    private boolean passFilter(Alignment alignment) {
        return alignment.isMapped() &&
                (includeDuplicates || !alignment.isDuplicate()) &&
                alignment.getMappingQuality() >= minMappingQuality &&
                !alignment.isVendorFailedRead();
    }

    /**
     * Parse and "count" the alignment file.  The main method.
     * <p/>
     * This method is not thread safe due to the use of the member variable "buffer".
     *
     * @throws IOException
     */
    public synchronized void parse() throws IOException {

        int tolerance = (int) (windowSize * (Math.floor(extFactor / windowSize) + 2));
        consumer.setSortTolerance(tolerance);

        AlignmentQueryReader reader = null;
        CloseableIterator<Alignment> iter = null;

        String lastChr = "";
        ReadCounter counter = null;

        WigWriter wigWriter = null;
        if (wigFile != null) {
            wigWriter = new WigWriter(wigFile, windowSize);
        }

        try {

            if (interval == null) {
                reader = getReader(alignmentFile, false);
                iter = reader.iterator();
            } else {
                reader = getReader(alignmentFile, true);
                iter = reader.query(interval.getChr(), interval.getStart(), interval.getEnd(), false);
            }

            while (iter != null && iter.hasNext()) {
                Alignment alignment = iter.next();
                if (passFilter(alignment)) {

                    totalCount++;

                    String alignmentChr = alignment.getChr();

                    // Close all counters with position < alignment.getStart()
                    if (alignmentChr.equals(lastChr)) {
                        if (counter != null) {
                            counter.closeBucketsBefore(alignment.getAlignmentStart() - tolerance, wigWriter);
                        }
                    } else {
                        if (counter != null) {
                            counter.closeBucketsBefore(Integer.MAX_VALUE, wigWriter);
                        }
                        counter = new ReadCounter(alignmentChr);
                        lastChr = alignmentChr;
                    }

                    AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
                    if (blocks != null) {
                        int lastBlockEnd = -1;
                        for (AlignmentBlock block : blocks) {

                            if (!block.isSoftClipped()) {

                                byte[] bases = block.getBases();
                                int blockStart = block.getStart();
                                int adjustedStart = block.getStart();
                                int adjustedEnd = block.getEnd();

                                if (alignment.isNegativeStrand()) {
                                    adjustedStart = Math.max(0, adjustedStart - extFactor);
                                } else {
                                    adjustedEnd += extFactor;
                                }

                                if (interval != null) {
                                    adjustedStart = Math.max(interval.getStart() - 1, adjustedStart);
                                    adjustedEnd = Math.min(interval.getEnd(), adjustedEnd);
                                }

                                for (int pos = adjustedStart; pos < adjustedEnd; pos++) {
                                    byte base = 0;
                                    int baseIdx = pos - blockStart;
                                    if (bases != null && baseIdx >= 0 && baseIdx < bases.length) {
                                        base = bases[baseIdx];
                                    }
                                    int idx = pos - blockStart;
                                    byte quality = (idx >= 0 && idx < block.qualities.length) ?
                                            block.qualities[pos - blockStart] : (byte) 0;
                                    counter.incrementCount(pos, base, quality);
                                }

                                lastBlockEnd = block.getEnd();
                            }
                        }
                    } else {
                        int adjustedStart = alignment.getAlignmentStart();
                        int adjustedEnd = alignment.getAlignmentEnd();

                        if (alignment.isNegativeStrand()) {
                            adjustedStart = Math.max(0, adjustedStart - extFactor);
                        } else {
                            adjustedEnd += extFactor;
                        }

                        if (interval != null) {
                            adjustedStart = Math.max(interval.getStart() - 1, adjustedStart);
                            adjustedEnd = Math.min(interval.getEnd(), adjustedEnd);
                        }


                        for (int pos = adjustedStart; pos < adjustedEnd; pos++) {
                            counter.incrementCount(pos, (byte) 0, (byte) 0);
                        }
                    }

                }

            }
            consumer.setAttribute("totalCount", String.valueOf(totalCount));
            consumer.parsingComplete();

        } catch (Exception e) {
            e.printStackTrace();
        } finally {

            if (counter != null) {
                counter.closeBucketsBefore(Integer.MAX_VALUE, wigWriter);
            }
            if (iter != null) {
                iter.close();
            }
            if (reader != null) {
                reader.close();
            }
            if (wigWriter != null) {
                wigWriter.close();
            }

        }
    }

    private AlignmentQueryReader getReader(String alignmentFile, boolean b) throws IOException {

        boolean isList = alignmentFile.indexOf(",") > 0;
        if (isList) {
            String[] tokens = alignmentFile.split(",");
            List<AlignmentQueryReader> readers = new ArrayList(tokens.length);
            for (String f : tokens) {
                readers.add(AlignmentReaderFactory.getReader(f, b));
            }
            return new MergedAlignmentReader(readers);
        } else {
            if (!FileUtils.isRemote(alignmentFile)) {
                File f = new File(alignmentFile);
                if (f.isDirectory()) {
                    List<AlignmentQueryReader> readers = new ArrayList();
                    for (File file : f.listFiles(new AlignmentFileFilter())) {
                        readers.add(AlignmentReaderFactory.getReader(file.getAbsolutePath(), b));
                    }
                    return new MergedAlignmentReader(readers);
                }
            }
            return AlignmentReaderFactory.getReader(alignmentFile, b);
        }

    }

    class ReadCounter {

        String chr;
        /**
         * Map of window index -> counter
         */
        TreeMap<Integer, Counter> counts = new TreeMap();

        ReadCounter(String chr) {
            this.chr = chr;
        }

        /**
         * @param position - genomic poistion
         * @param base     - nucleotide
         * @param quality  - base quality of call
         */
        void incrementCount(int position, byte base, byte quality) {
            final Counter counter = getCounterForPosition(position);
            counter.increment(position, base, quality);
        }


        private Counter getCounterForPosition(int position) {
            int idx = position / windowSize;
            return getCounter(idx);
        }

        private Counter getCounter(int idx) {
            if (!counts.containsKey(idx)) {
                int counterStartPosition = idx * windowSize;
                int counterEndPosition = counterStartPosition + windowSize;
                counts.put(idx, new Counter(chr, counterStartPosition, counterEndPosition));
            }
            final Counter counter = counts.get(idx);
            return counter;
        }


        /**
         * Close (finalize) all buckets before the given position.  Called when we are sure this position will not be
         * visited again.
         *
         * @param position - genomic position
         */
        void closeBucketsBefore(int position, WigWriter wigWriter) {
            List<Integer> bucketsToClose = new ArrayList();

            Integer bucket = position / windowSize;
            for (Map.Entry<Integer, Counter> entry : counts.entrySet()) {
                if (entry.getKey() < bucket) {

                    final Counter counter = entry.getValue();
                    int totalCount = counter.getCount();

                    // Divide total count by window size.  This is the average count per
                    // base over the window,  so for example 30x coverage remains 30x irrespective of window size.
                    int bucketStartPosition = entry.getKey() * windowSize;
                    int bucketEndPosition = bucketStartPosition + windowSize;
                    if (genome != null) {
                        Chromosome chromosome = genome.getChromosome(chr);
                        if (chromosome != null) {
                            bucketEndPosition = Math.min(bucketEndPosition, chromosome.getLength());
                        }
                    }
                    int bucketSize = bucketEndPosition - bucketStartPosition;

                    buffer[0] = ((float) totalCount) / bucketSize;

                    if (strandOption > 0) {
                        buffer[1] = ((float) counter.getCount()) / bucketSize;
                    }

                    consumer.addData(chr, bucketStartPosition, bucketEndPosition, buffer, null);

                    if (wigWriter != null) {
                        wigWriter.addData(chr, bucketStartPosition, bucketEndPosition, buffer);
                    }


                    bucketsToClose.add(entry.getKey());
                }
            }

            for (Integer key : bucketsToClose) {
                counts.remove(key);
            }


        }

    }


    /**
     * Events
     * base mismatch
     * translocation
     * insertion (small)
     * insertion (large, rearrangment)
     * deletion  (small)
     * deletion  (large, rearrangment)
     */
    class Counter {
        int count = 0;
        int negCount = 0;
        int qualityCount = 0;

        float pairedCount = 0;
        float mismatchCount = 0;


        String chr;
        int start;
        int end;
        byte[] ref;
        int[] baseCount;


        Counter(String chr, int start, int end) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            baseCount = new int[end - start];
        }


        int getCount() {
            return count;
        }

        int getNegCount() {
            return negCount;
        }


        void incrementNeg() {
            negCount++;
        }


        // frac should be between 0 znd 1
        void incrementPairedCount(float frac) {
            pairedCount += frac;
        }


        public int[] getBaseCount() {
            return baseCount;
        }


        void increment(int position, byte base, byte quality) {

            // Qualities of 2 or less => no idea what this base is
            //if (quality <= 2) {
            //    return;
            //}

            int offset = position - start;
            baseCount[offset]++;

            if (ref != null && ref.length > offset) {
                byte refBase = ref[offset];
                if (refBase != base) {
                    mismatchCount += quality;
                }
            }
            count++;
            qualityCount += quality;
        }


    }


    /**
     * Creates a vary step wig file
     */
    class WigWriter {


        String lastChr = null;
        int lastPosition = 0;
        int step;
        int span;
        PrintWriter pw;

        WigWriter(File file, int step) throws IOException {
            this.step = step;
            this.span = step;
            pw = new PrintWriter(new FileWriter(file));
        }

        public void addData(String chr, int start, int end, float[] data) {

            for (int i = 0; i < data.length; i++) {
                if (Float.isNaN(data[i])) {
                    return;
                }
            }

            if (genome.getChromosome(chr) == null) {
                return;
            }
            if (end <= start) {     // Not sure why or how this could happen
                return;
            }

            int dataSpan = end - start;

            if (chr == null || !chr.equals(lastChr) || dataSpan != span) {
                span = dataSpan;
                outputStepLine(chr, start + 1);
            }

            pw.print(start + 1);
            for (int i = 0; i < data.length; i++) {
                System.out.print("\t" + data[i]);
            }
            pw.println();

            lastPosition = start;
            lastChr = chr;

        }

        private void close() {
            pw.close();

        }

        private void outputStepLine(String chr, int start) {
            pw.println("variableStep chrom=" + chr + " span=" + span);
        }

    }


}
