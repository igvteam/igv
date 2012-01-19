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
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.MergedAlignmentReader;
import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.ui.filefilters.AlignmentFileFilter;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.Utilities;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

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
     * 0x01 - Strand agnostic (both)
     * 0x02 - Separate strands by first read/second read,
     * instead of by positive/negative (if this flag is unset)
     * 0x04 - Positive/First strand
     * 0x08 - Negative/Second strand
     * TODO Make these static constants or enum
     */
    private int strandOption = 0x01;
    private static final int NUM_STRANDS = 2;
    private static final int STRAND_TOTAL = 0x01;
    private static final int FIRST_IN_PAIR = 0x02;
    private static final int STRAND_ZERO = 0x04;
    private static final int STRAND_ONE = 0x08;

    private static boolean outputZero;
    private static boolean outputOne;
    private static boolean outputTotal;
    private static boolean firstInPair;

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


    private final static Set<Byte> nucleotidesKeep = new HashSet<Byte>();

    static {
        byte[] tokeep = new byte[]{'A', 'C', 'G', 'T', 'N'};
        for (byte b : tokeep) {
            nucleotidesKeep.add(b);
        }
    }


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

        if (options != null) {
            parseOptions(options);
        }

        buffer = new float[Utilities.countFlags(strandOption)];

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
            } else if (opt.startsWith("sc=")) {
                String[] tmp = opt.split("=");
                strandOption = Integer.decode(tmp[1]);
            } else {
                System.out.println("Unknown coverage option: " + opt);
            }
        }
        outputZero = (strandOption & STRAND_ZERO) > 0;
        outputOne = (strandOption & STRAND_ONE) > 0;
        outputTotal = (strandOption & STRAND_TOTAL) > 0;
        firstInPair = (strandOption & FIRST_IN_PAIR) > 0;
    }


    // TODO -- command-line options to override all of these checks
    private boolean passFilter(Alignment alignment) {
        boolean pairingInfo = !firstInPair || (alignment.getFirstOfPairStrand() != Strand.NONE);

        return alignment.isMapped() && pairingInfo &&
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
                    //Sort into the read strand or first-in-pair strand,
                    //depending on input flag
                    Strand strand;
                    if (firstInPair) {
                        strand = alignment.getFirstOfPairStrand();
                    } else {
                        strand = alignment.getReadStrand();
                    }
                    if (strand.equals(Strand.NONE)) {
                        //TODO move this into passFilter, or move passFilter here
                        continue;
                    }
                    boolean readNegStrand = alignment.isNegativeStrand();

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

                                if (readNegStrand) {
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
                                    counter.incrementCount(pos, base, quality, strand);
                                }

                                lastBlockEnd = block.getEnd();
                            }
                        }
                    } else {
                        int adjustedStart = alignment.getAlignmentStart();
                        int adjustedEnd = alignment.getAlignmentEnd();

                        if (readNegStrand) {
                            adjustedStart = Math.max(0, adjustedStart - extFactor);
                        } else {
                            adjustedEnd += extFactor;
                        }

                        if (interval != null) {
                            adjustedStart = Math.max(interval.getStart() - 1, adjustedStart);
                            adjustedEnd = Math.min(interval.getEnd(), adjustedEnd);
                        }


                        for (int pos = adjustedStart; pos < adjustedEnd; pos++) {
                            counter.incrementCount(pos, (byte) 0, (byte) 0, strand);
                        }
                    }

                }

            }
            consumer.setAttribute("totalCount", String.valueOf(totalCount));

            //TODO Mostly just for testing
//            int[] strandCounts = new int[NUM_STRANDS];
//            for(Map.Entry<Integer, Counter> counts: counter.counts.entrySet()){
//                int[][] baseCount = counts.getValue().getBaseCount();
//                for(int ii=0; ii < strandCounts.length; ii++){
//                    int[] tmp = baseCount[ii];
//                    if( tmp == null) continue;
//                    for(int jj=0; jj < tmp.length; jj++){
//                        strandCounts[ii] += tmp[jj];
//                    }
//                }
//            }
            //consumer.setAttribute("positiveStrandCount", "" + strandCounts[Strand.POSITIVE.ordinal()]);
            //consumer.setAttribute("negativeStrandCount", "" + strandCounts[Strand.NEGATIVE.ordinal()]);
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
         * @param position - genomic position
         * @param base     - nucleotide
         * @param quality  - base quality of call
         * @param strand   - which strand to increment count. Should be POSITIVE or NEGATIVE
         */
        void incrementCount(int position, byte base, byte quality, Strand strand) {
            final Counter counter = getCounterForPosition(position);
            int strandNum = strand.equals(Strand.POSITIVE) ? 0 : 1;
            counter.increment(position, base, quality, strandNum);
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

                    final Counter counter = entry.getValue();

                    boolean[] outCols = new boolean[]{outputZero, outputOne};
                    int col = 0;

                    //Output aggregated information
                    if (outputTotal) {
                        buffer[col] = ((float) counter.getTotalCounts()) / bucketSize;
                        col++;
                    }
                    //Output strand specific information, if applicable
                    for (int ii = 0; ii < NUM_STRANDS; ii++) {
                        if (outCols[ii]) {
                            buffer[col] = ((float) counter.getCount(ii)) / bucketSize;
                            col++;
                        }

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
     * Class for counting nucleotides and strands over an interval.
     */

    class Counter {
        /**
         * The total number of counts on this Counter. Will
         * be the sum of baseCount over the second index.
         */
        int[] strandCount = new int[NUM_STRANDS];
        int totalCount = 0;
        int qualityCount = 0;

        String chr;
        int start;
        int end;
        //byte[] ref;
        /**
         * The number of times a particular base has been encountered (ie # of reads of that base)
         */
        //int[][] baseCount = new int[NUM_STRANDS][];

        private Map<Byte, Integer> baseTypeCounts = new HashMap();

        Counter(String chr, int start, int end) {
            this.chr = chr;
            this.start = start;
            this.end = end;
        }


        int getCount(int strand) {
            return strandCount[strand];
        }

        // TODO -- do we need to expose this implementation detail?
//        public int[][] getBaseCount() {
//            return baseCount;
//        }


        void increment(int position, byte base, byte quality, int strand) {
            //int offset = position - start;
            //Lazy initialization
//            if (baseCount[strand] == null) {
//                baseCount[strand] = new int[this.end - this.start];
//            }
            //baseCount[strand][offset]++;

            incrementNucleotide(base);

            strandCount[strand]++;
            totalCount++;
            qualityCount += quality;
        }

        /**
         * Increment the nucleotide counts.
         *
         * @param base 65, 67, 71, 84, 78
         *             aka A, C, G, T, N (upper case).
         *             Anything else is stored as 0
         */
        private void incrementNucleotide(byte base) {
            if (!nucleotidesKeep.contains(base)) {
                base = 0;
            }

            int orig = 0;
            if (baseTypeCounts.containsKey(base)) {
                orig = baseTypeCounts.get(base);
            }
            baseTypeCounts.put(base, orig + 1);
        }


        public int getTotalCounts() {
            return this.totalCount;
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
                pw.print("\t" + data[i]);
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
