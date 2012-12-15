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
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.tools.parsers.DataConsumer;

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

    /*
     * Output data from each strand separately (as opposed to combining them)
     * using the read strand.
     */
    static final int STRANDS_BY_READ = 0x01;

    /*
     * Output strand data separately by first-in-pair
     */
    static final int STRANDS_BY_FIRST_IN_PAIR = 0x02;

    /**
     * Output strand data separately by second-in-pair
     */
    static final int STRANDS_BY_SECOND_IN_PAIR = 0x04;

    /**
     * Output counts of each base. Whether the data will be output
     * for each strand separately is determined by STRAND_SEPARATE
     * by
     */
    static final int BASES = 0x08;

    public static final int INCLUDE_DUPS = 0x20;
    public static final int PAIRED_COVERAGE = 0x40;
    private boolean outputSeparate;
    private boolean firstInPair;
    private boolean secondInPair;
    private boolean outputBases;

    private static final int[] output_strands = new int[]{0, 1};

    public static final int NUM_STRANDS = output_strands.length;

    /**
     * Extension factor.  Reads are extended by this amount from the 3' end before counting.   The purpose is to yield
     * an approximate count of fragment "coverage", as opposed to read coverage.  If used, the value should be set to
     * extFactor = averageFragmentSize - averageReadLength
     */
    private int extFactor;

    /**
     * 5' "pre" extension factor.  Read is extended by this amount from the 5' end of the read
     */
    private int preExtFactor;

    /**
     * 5' "post" extension factor.  Essentially, replace actual read length by this amount.
     */
    private int postExtFactor;

    /**
     * Flag to control treatment of duplicates.  If true duplicates are counted.  The default value is false.
     */
    private boolean includeDuplicates = false;

    /**
     * If true, coverage is computed based from properly paired reads, counting entire insert.
     */
    private boolean pairedCoverage = false;

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
    private Locus queryInterval;

    /**
     * Data buffer to pass data to the "consumer" (preprocessor).
     */
    private float[] buffer;


    private boolean computeTDF = true;


    private final static Set<Byte> nucleotidesKeep = new HashSet<Byte>();
    private final static byte[] nucleotides = new byte[]{'A', 'C', 'G', 'T', 'N'};

    static {
        for (byte b : nucleotides) {
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
     * @param queryString   - Locus query string, such as 1:1-1000. Only count the queried region. Set to null for entire genome
     * @param minMapQual    - Minimum mapping quality to include
     * @param countFlags    - Combination of flags for BASES, STRAND_SEPARATE, INCLUDE_DUPES, FIRST_IN_PAIR
     */
    public CoverageCounter(String alignmentFile,
                           DataConsumer consumer,
                           int windowSize,
                           int extFactor,
                           File wigFile,
                           Genome genome,
                           String queryString,
                           int minMapQual,
                           int countFlags) {
        this.alignmentFile = alignmentFile;
        this.consumer = consumer;
        this.windowSize = windowSize;
        this.extFactor = extFactor;
        this.wigFile = wigFile;
        this.genome = genome;

        parseOptions(queryString, minMapQual, countFlags);

        //Count the number of output columns. 1 or 2 if not outputting bases
        //5 or 10 if are.
        int multiplier = outputBases ? 5 : 1;
        int datacols = (outputSeparate ? 2 : 1) * multiplier;

        buffer = new float[datacols];
    }

    public void setPreExtFactor(int preExtFactor) {
        this.preExtFactor = preExtFactor;
    }

    public void setPosExtFactor(int postExtFactor) {
        this.postExtFactor = postExtFactor;
    }

    /**
     * Take additional optional command line arguments and parse them
     *
     * @param queryString
     * @param minMapQual
     * @param countFlags
     */
    private void parseOptions(String queryString, int minMapQual, int countFlags) {
        if (queryString != null) {
            this.queryInterval = new Locus(queryString);
        }
        this.minMappingQuality = minMapQual;
        outputSeparate = (countFlags & STRANDS_BY_READ) > 0;
        firstInPair = (countFlags & STRANDS_BY_FIRST_IN_PAIR) > 0;
        secondInPair = (countFlags & STRANDS_BY_SECOND_IN_PAIR) > 0;
        outputSeparate |= firstInPair || secondInPair;
        if (firstInPair && secondInPair) {
            throw new IllegalArgumentException("Can't set both first and second in pair");
        }
        outputBases = (countFlags & BASES) > 0;
        includeDuplicates = (countFlags & INCLUDE_DUPS) > 0;
        pairedCoverage = (countFlags & PAIRED_COVERAGE) > 0;
    }


    // TODO -- command-line options to override all of these checks
    private boolean passFilter(Alignment alignment) {

        // If the first-in-pair or second-in-pair option is selected test that we have that information, otherwise
        // alignment is filtered.
        boolean pairingInfo = (!firstInPair && !secondInPair) ||
                (alignment.getFirstOfPairStrand() != Strand.NONE);

        // For paired coverage, see if the alignment is properly paired, and if it is the "leftmost" alignment
        // (to prevent double-counting the pair).
        if (pairedCoverage) {
            ReadMate mate = alignment.getMate();
            if (!alignment.isProperPair() || alignment.getMate() == null || alignment.getStart() > mate.getStart()) {
                return false;
            }
            if (Math.abs(alignment.getInferredInsertSize()) > 10000) {
                System.out.println("Very large insert size: " + Math.abs(alignment.getInferredInsertSize()) +
                        " for read " + alignment.getReadName() + ".  Skipped.");
                return false;
            }
        }

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

        int maxExtFactor = Math.max(extFactor, Math.max(preExtFactor, postExtFactor));

        int tolerance = (int) (windowSize * (Math.floor(maxExtFactor / windowSize) + 2));
        consumer.setSortTolerance(tolerance);

        AlignmentReader reader = null;
        CloseableIterator<Alignment> iter = null;

        String lastChr = "";
        ReadCounter counter = null;

        WigWriter wigWriter = null;
        if (wigFile != null) {
            wigWriter = new WigWriter(wigFile, windowSize);
        }

        try {

            if (queryInterval == null) {
                reader = AlignmentReaderFactory.getReader(alignmentFile, false);
                iter = reader.iterator();
            } else {
                reader = AlignmentReaderFactory.getReader(alignmentFile, true);
                iter = reader.query(queryInterval.getChr(), queryInterval.getStart() - 1, queryInterval.getEnd(), false);
            }

            while (iter != null && iter.hasNext()) {
                Alignment alignment = iter.next();
                if (passFilter(alignment)) {
                    //Sort into the read strand or first-in-pair strand,
                    //depending on input flag. Note that this can
                    //be very unreliable depending on data
                    Strand strand;
                    if (firstInPair) {
                        strand = alignment.getFirstOfPairStrand();
                    } else if (secondInPair) {
                        strand = alignment.getSecondOfPairStrand();
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
                    } else {  // New chromosome
                        if (counter != null) {
                            counter.closeBucketsBefore(Integer.MAX_VALUE, wigWriter);
                        }
                        counter = new ReadCounter(alignmentChr);
                        lastChr = alignmentChr;
                    }

                    AlignmentBlock[] blocks = alignment.getAlignmentBlocks();

                    if (blocks != null && !pairedCoverage) {
                        int lastBlockEnd = -1;
                        for (AlignmentBlock block : blocks) {

                            if (!block.isSoftClipped()) {

                                byte[] bases = block.getBases();
                                int blockStart = block.getStart();
                                int blockEnd = block.getEnd();


                                int adjustedStart = block.getStart();
                                int adjustedEnd = block.getEnd();


                                if (preExtFactor > 0) {
                                    if (readNegStrand) {
                                        adjustedEnd = blockEnd + preExtFactor;
                                    } else {
                                        adjustedStart = Math.max(0, blockStart - preExtFactor);
                                    }
                                }

                                // If both postExtFactor and extFactor are specified, postExtFactor takes precedence
                                //
                                if (postExtFactor > 0) {
                                    if (readNegStrand) {
                                        adjustedStart = Math.max(0, blockEnd - postExtFactor);
                                    } else {
                                        adjustedEnd = blockStart + postExtFactor;
                                    }

                                } else if (extFactor > 0) {
                                    // Standard extension option -- extend read on 3' end
                                    if (readNegStrand) {
                                        adjustedStart = Math.max(0, adjustedStart - extFactor);
                                    } else {
                                        adjustedEnd += extFactor;
                                    }
                                }


                                if (queryInterval != null) {
                                    adjustedStart = Math.max(queryInterval.getStart() - 1, adjustedStart);
                                    adjustedEnd = Math.min(queryInterval.getEnd(), adjustedEnd);
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
                        int adjustedEnd = pairedCoverage ?
                                adjustedStart + Math.abs(alignment.getInferredInsertSize()) :
                                alignment.getAlignmentEnd();

                        if (readNegStrand) {
                            adjustedStart = Math.max(0, adjustedStart - extFactor);
                        } else {
                            adjustedEnd += extFactor;
                        }

                        if (queryInterval != null) {
                            adjustedStart = Math.max(queryInterval.getStart() - 1, adjustedStart);
                            adjustedEnd = Math.min(queryInterval.getEnd(), adjustedEnd);
                        }


                        for (int pos = adjustedStart; pos < adjustedEnd; pos++) {
                            counter.incrementCount(pos, (byte) 'N', (byte) 0, strand);
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


    /**
     * The names of tracks which will be created by this parser
     *
     * @param prefix String to be prepended to each track name
     * @return
     */
    public String[] getTrackNames(String prefix) {
        if (prefix == null) {
            prefix = "";
        }
        String[] trackNames = new String[this.buffer.length];
        String[] strandArr;
        if (outputSeparate) {
            strandArr = new String[]{"Positive Strand", "Negative Strand"};
        } else {
            strandArr = new String[]{"Combined Strands"};
        }
        int col = 0;
        for (String sA : strandArr) {
            if (outputBases) {
                for (Byte n : nucleotides) {
                    trackNames[col] = prefix + " " + sA + " " + new String(new byte[]{n});
                    col++;
                }
            } else {
                trackNames[col] = prefix + " " + sA;
                col++;
            }
        }
        return trackNames;
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
            counter.increment(base, quality, strandNum);
        }


        private Counter getCounterForPosition(int position) {
            int idx = position / windowSize;
            return getCounter(idx);
        }

        private Counter getCounter(int idx) {
            if (!counts.containsKey(idx)) {
                counts.put(idx, new Counter());
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

            int bucket = position / windowSize;
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

                    int col = 0;

                    //Not outputting base info, just totals
                    if (!outputBases) {
                        if (outputSeparate) {
                            //Output strand specific information, if applicable
                            for (int strandNum : output_strands) {
                                buffer[col] = ((float) counter.getCount(strandNum)) / bucketSize;
                                col++;
                            }

                        } else {
                            buffer[col] = ((float) counter.getTotalCounts()) / bucketSize;
                            col++;
                        }

                        //Output counts of each base
                    } else {
                        if (outputSeparate) {
                            for (int strandNum : output_strands) {
                                for (byte base : nucleotides) {
                                    buffer[col] = ((float) counter.getBaseCount(base, strandNum)) / bucketSize;
                                    col++;
                                }
                            }
                        } else {
                            for (byte base : nucleotides) {
                                buffer[col] = ((float) counter.getBaseCount(base)) / bucketSize;
                                col++;
                            }
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
        int[] strandCount;
        int totalCount = 0;
        //int qualityCount = 0;

        //String chr;
        //int start;
        //int end;
        //byte[] ref;
        /**
         * The number of times a particular base has been encountered (ie # of reads of that base)
         */
        //int[][] baseCount = new int[NUM_STRANDS][];

        private Map<Byte, Integer>[] baseTypeCounts;

        Counter() {
            //this.chr = chr;
            //this.start = start;
            //this.end = end;

            if (outputBases) {
                baseTypeCounts = new HashMap[NUM_STRANDS];
                for (int ii = 0; ii < NUM_STRANDS; ii++) {
                    baseTypeCounts[ii] = new HashMap<Byte, Integer>();
                }
            }

            if (outputSeparate) {
                strandCount = new int[NUM_STRANDS];
            }
        }


        int getCount(int strand) {
            return strandCount[strand];
        }

        public int getTotalCounts() {
            return totalCount;
        }

        void increment(byte base, byte quality, int strand) {

            if (outputBases) {
                incrementNucleotide(base, strand);
            }

            if (outputSeparate) {
                this.strandCount[strand]++;
            }

            this.totalCount++;
            //this.qualityCount += quality;
        }

        /**
         * Increment the nucleotide counts.
         *
         * @param base   65, 67, 71, 84, 78
         *               aka A, C, G, T, N (upper case).
         *               Anything else is stored as 0
         * @param strand index of strand, 0 for positive and 1 for negative
         */
        private void incrementNucleotide(byte base, int strand) {
            Map<Byte, Integer> btc = baseTypeCounts[strand];
            if (!nucleotidesKeep.contains(base)) {
                base = 0;
            }

            int orig = 0;
            if (btc.containsKey(base)) {
                orig = btc.get(base);
            }
            btc.put(base, orig + 1);
        }

        public int getBaseCount(byte base, int strand) {
            return baseTypeCounts[strand].containsKey(base) ? baseTypeCounts[strand].get(base) : 0;
        }

        public int getBaseCount(byte base) {
            int count = 0;
            for (int strand = 0; strand < NUM_STRANDS; strand++) {
                count += getBaseCount(base, strand);
            }
            return count;
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
                outputHeader(chr, start + 1);
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

        private void outputHeader(String chr, int start) {
            //Write label column
            String labels = "Pos";
            for (String s : getTrackNames("")) {
                labels += "," + s;
            }
            pw.println("#Columns: " + labels);
            pw.println("variableStep chrom=" + chr + " span=" + span);
        }

    }


}
