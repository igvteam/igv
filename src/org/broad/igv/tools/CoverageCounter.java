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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tools;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.SequenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.Strand;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.SamQueryReaderFactory;
import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.util.stats.Histogram;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;

/**
 *   TODO -- normalize option
 *
 -n           Normalize the count by the total number of reads. This option
 multiplies each count by (1,000,000 / total # of reads). It
 is useful when comparing multiple chip-seq experiments when
 the absolute coverage depth is not important.
 */

/**
 *
 */
public class CoverageCounter {

    private String alignmentFile;
    private File tdfFile;
    private DataConsumer consumer;
    private float[] buffer;
    private int windowSize = 1;
    // TODO -- make mapping qulaity a parameter
    private int minMappingQuality = 0;
    private int strandOption = -1;
    private int extFactor;
    private int totalCount = 0;
    private File wigFile = null;
    private WigWriter wigWriter = null;
    private Genome genome;

    WigWriter mismatchWigWriter;
    WigWriter isizeWigWriter;
    WigWriter orientationWigWriter;
    WigWriter indelWigWriter;
    WigWriter nomateWigWriter;

    //private float meanInsertSize;
    //private float stdDevInsertSize;
    private boolean computeTDF = true;

    private Histogram coverageHistogram;
    private static final double LOG_1__1 = 0.09531018;


    public CoverageCounter(String alignmentFile,
                           DataConsumer consumer,
                           int windowSize,
                           int extFactor,
                           File tdfFile,  // For reference
                           File wigFile,
                           Genome genome,
                           int strandOption,
                           String options) {
        this.alignmentFile = alignmentFile;
        this.tdfFile = tdfFile;
        this.consumer = consumer;
        this.windowSize = windowSize;
        this.extFactor = extFactor;
        this.wigFile = wigFile;
        this.genome = genome;
        this.strandOption = strandOption;
        buffer = strandOption < 0 ? new float[1] : new float[2];

        if (options != null) {
            parseOptions(options);
        }
    }

    private void parseOptions(String options) {
        String[] opts = options.split(",");
        for (String opt : opts) {

            if (opt.equals("i")) {
                computeISize();
            } else if (opt.equals("o")) {
                computeOrientation();
            } else if (opt.equals("m")) {
                computeMismatch();
            } else if (opt.equals("d")) {
                computeIndel();
            } else if (opt.equals("n")) {
                computeNomate();
            } else if (opt.equals("h")) {
                coverageHistogram = new Histogram(1000);
            } else {
                System.out.println("Unknown coverage option: " + opt);
            }
        }
    }

    public void computeMismatch() {
        try {
            mismatchWigWriter = new WigWriter(new File(tdfFile.getAbsolutePath() + ".mismatch.wig"), windowSize, false);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    public void computeISize() { //float meanInsertSize, float stdDevInsertSize) {
        try {
            //this.meanInsertSize = meanInsertSize;
            //this.stdDevInsertSize = stdDevInsertSize;
            isizeWigWriter = new WigWriter(new File(tdfFile.getAbsolutePath() + ".isize.wig"), windowSize, false);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    public void computeOrientation() {
        try {
            orientationWigWriter = new WigWriter(new File(tdfFile.getAbsolutePath() + ".orientation.wig"), windowSize, false);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    public void computeIndel() {
        try {
            indelWigWriter = new WigWriter(new File(tdfFile.getAbsolutePath() + ".indel.wig"), windowSize, false);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    public void computeNomate() {
        try {
            nomateWigWriter = new WigWriter(new File(tdfFile.getAbsolutePath() + ".nomate.wig"), windowSize, false);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


    // TODO -- options to ovveride all of these checks

    private boolean passFilter(Alignment alignment) {

        if (strandOption > 0 && alignment.getFragmentStrand(strandOption) == Strand.NONE) {
            return false;
        }
        return alignment.isMapped() &&
                !alignment.isDuplicate() &&
                alignment.getMappingQuality() >= minMappingQuality &&
                !alignment.isVendorFailedRead() &&
                !alignment.isDuplicate();
    }

    public void parse() {

        int tolerance = (int) (windowSize * (Math.floor(extFactor / windowSize) + 2));
        consumer.setSortTolerance(tolerance);

        AlignmentQueryReader reader = null;
        CloseableIterator<Alignment> iter = null;

        String lastChr = "";
        ReadCounter counter = null;


        try {

            if (wigFile != null) {
                wigWriter = new WigWriter(wigFile, windowSize, false);
            }


            reader = SamQueryReaderFactory.getReader(alignmentFile, false);
            iter = reader.iterator();

            while (iter != null && iter.hasNext()) {
                Alignment alignment = iter.next();
                if (passFilter(alignment)) {

                    totalCount++;

                    String alignmentChr = alignment.getChr();

                    // Close all counters with position < alignment.getStart()
                    if (alignmentChr.equals(lastChr)) {
                        if (counter != null) {
                            counter.closeBucketsBefore(alignment.getAlignmentStart() - tolerance);
                        }
                    } else {
                        if (counter != null) {
                            counter.closeBucketsBefore(Integer.MAX_VALUE);
                        }
                        counter = new ReadCounter(alignmentChr);
                        lastChr = alignmentChr;
                    }

                    if (alignment.isPaired()) {

                        counter.incrementPairedCount (alignment.getStart());

                        ReadMate mate = alignment.getMate();
                        boolean mateMapped = mate != null && mate.isMapped();
                        boolean sameChromosome = mateMapped && mate.getChr().equals(alignment.getChr());

                        if (isizeWigWriter != null && mateMapped && sameChromosome) {
                            int isize = Math.abs(alignment.getInferredInsertSize());
                            float score = Math.min(900, isize);
                            if (isize > 900) {
                                score += (float) (Math.log(isize - 900) / LOG_1__1);
                            }
                            // Plot looks better if we record this for a point (start pos) rather than the alignment extent
                            isizeWigWriter.addData(alignment.getChr(), alignment.getStart(), alignment.getStart() + 1, score);
                        }

                        //if(interChrWigWriter != null && mateMapped && !sameChromosome) {
                        //
                        //}

                        if (orientationWigWriter != null && mateMapped && sameChromosome) {
                            String oStr = alignment.getPairOrientation();
                            int orientation = getOrientation(oStr);
                            counter.incrementOrientation(alignment.getStart(), orientation);
                        }

                        if (nomateWigWriter != null && !mateMapped) {
                            counter.incrementNomate(alignment.getStart());
                        }
                    }


                    AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
                    if (blocks != null) {
                        int lastBlockEnd = -1;
                        for (AlignmentBlock block : blocks) {

                            if (indelWigWriter != null && lastBlockEnd >= 0) {
                                int s = block.getStart();
                                counter.incrementIndel(lastBlockEnd, s);
                            }

                            byte[] bases = block.getBases();
                            int blockStart = block.getStart();
                            int adjustedStart = block.getStart();
                            int adjustedEnd = block.getEnd();
                            if (alignment.isNegativeStrand()) {
                                adjustedStart = Math.max(0, adjustedStart - extFactor);
                            } else {
                                adjustedEnd += extFactor;
                            }

                            for (int pos = adjustedStart; pos < adjustedEnd; pos++) {
                                byte base = 0;
                                int baseIdx = pos - blockStart;
                                if (bases != null && baseIdx >= 0 && baseIdx < bases.length) {
                                    base = bases[baseIdx];
                                }
                                counter.incrementCount(pos, base);
                            }

                            lastBlockEnd = block.getEnd();
                        }
                    } else {
                        int adjustedStart = alignment.getAlignmentStart();
                        int adjustedEnd = alignment.getAlignmentEnd();
                        if (alignment.isNegativeStrand()) {
                            adjustedStart = Math.max(0, adjustedStart - extFactor);
                        } else {
                            adjustedEnd += extFactor;
                        }

                        for (int pos = adjustedStart; pos < adjustedEnd; pos++) {
                            counter.incrementCount(pos, (byte) 0);
                        }
                    }

                    if (indelWigWriter != null) {
                        for (AlignmentBlock block : alignment.getInsertions()) {
                            counter.incrementIndel(block.getStart() - 1, block.getStart());
                        }

                    }
                }

            }
        }

        catch (
                Exception e
                )

        {
            e.printStackTrace();
        }

        finally

        {
            if (counter != null) {
                counter.closeBucketsBefore(Integer.MAX_VALUE);
            }

            consumer.setAttribute("totalCount", String.valueOf(totalCount));
            consumer.parsingComplete();

            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

            if (iter != null) {
                iter.close();
            }

            if (wigWriter != null) {
                wigWriter.close();
            }

            if (mismatchWigWriter != null)
                mismatchWigWriter.close();

            if (isizeWigWriter != null) {
                isizeWigWriter.close();
            }
            if (orientationWigWriter != null) {
                orientationWigWriter.close();
            }
            if (indelWigWriter != null) {
                indelWigWriter.close();
            }
            if (nomateWigWriter != null) {
                nomateWigWriter.close();
            }
            if (coverageHistogram != null) {
                try {
                    PrintWriter pw = new PrintWriter(new FileWriter(tdfFile.getAbsolutePath() + ".hist.txt"));
                    coverageHistogram.print(pw);
                    pw.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
    }

    class ReadCounter {

        String chr;
        TreeMap<Integer, Counter> counts = new TreeMap();

        ReadCounter(String chr) {
            this.chr = chr;
        }

        void incrementCount(int position, byte base) {
            final Counter counter = getBucket(position);
            counter.increment(position, base);
        }


        void incrementIndel(int startPosition, int endPosition) {
            for (int position = startPosition; position < endPosition; position++) {
                final Counter counter = getBucket(position);
                counter.incrementIndel(position);
            }
        }

        void incrementNomate(int position) {
            final Counter counter = getBucket(position);
            counter.incrementNomate();
        }

        void incrementOrientation(int position, int orientation) {
            final Counter counter = getBucket(position);
            counter.incrementOrientation(orientation);
        }

        void incrementPairedCount(int position) {
            final Counter counter = getBucket(position);
            counter.incrementPairedCount();
        }



        public void incrementISize(Alignment alignment) {
            int startBucket = alignment.getStart() / windowSize;
            int endBucket = alignment.getAlignmentEnd() / windowSize;
            for (int bucket = startBucket; bucket <= endBucket; bucket++) {
                int bucketStartPosition = bucket * windowSize;
                int bucketEndPosition = bucketStartPosition + windowSize;
                if (!counts.containsKey(bucket)) {
                    counts.put(bucket, new Counter(chr, bucketStartPosition, bucketEndPosition));
                }
                final Counter counter = counts.get(bucket);
                counter.incrementISize(alignment.getInferredInsertSize());
            }
        }


        private Counter getBucket(int position) {
            int bucket = position / windowSize;
            int bucketStartPosition = bucket * windowSize;
            int bucketEndPosition = bucketStartPosition + windowSize;
            if (!counts.containsKey(bucket)) {
                counts.put(bucket, new Counter(chr, bucketStartPosition, bucketEndPosition));
            }
            final Counter counter = counts.get(bucket);
            return counter;
        }


        //void incrementNegCount(int position) {
        //    Integer bucket = position / windowSize;
        //    if (!counts.containsKey(bucket)) {
        //        counts.put(bucket, new Counter());
        //    }
        //    counts.get(bucket).incrementNeg();
        //}

        void closeBucketsBefore(int position) {
            List<Integer> bucketsToClose = new ArrayList();

            Integer bucket = position / windowSize;
            for (Map.Entry<Integer, Counter> entry : counts.entrySet()) {
                if (entry.getKey() < bucket) {

                    // Divide total count by window size.  This is the average count per
                    // base over the window,  so 30x coverage remains 30x irrespective of window size.
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
                    buffer[0] = ((float) counter.getCount()) / bucketSize;

                    if (strandOption > 0) {
                        buffer[1] = ((float) counter.getCount()) / bucketSize;
                    }

                    consumer.addData(chr, bucketStartPosition, bucketEndPosition, buffer, null);

                    if (mismatchWigWriter != null) {
                        float mismatch = counter.getMaxMismatchFraction();
                        mismatchWigWriter.addData(chr, bucketStartPosition, bucketEndPosition, mismatch);
                    }

                    if (indelWigWriter != null) {
                        float indelCount = counter.getMaxIndelFraction();
                        indelWigWriter.addData(chr, bucketStartPosition, bucketEndPosition, indelCount);
                    }

                    if (nomateWigWriter != null) {
                        float nomateCount = counter.getNomateFraction();
                        nomateWigWriter.addData(chr, bucketStartPosition, bucketEndPosition, nomateCount);
                    }
                    if (orientationWigWriter != null) {
                        float posOrientation = counter.getPosOrientationFraction();
                        if (posOrientation > 0) {
                            orientationWigWriter.addData(chr, bucketStartPosition, bucketEndPosition, posOrientation);
                        }
                        float negOrientation = counter.getNegOrientationFraction();
                        if (negOrientation > 0) {
                            orientationWigWriter.addData(chr, bucketStartPosition, bucketEndPosition, negOrientation);
                        }

                    }

                    //if (isizeWigWriter != null) {
                    //float isizeScore = counter.getISizeFraction();
                    //float isizeScore = counter.getAvgIsize();
                    //isizeWigWriter.addData(chr, bucketStartPosition, bucketEndPosition, isizeScore);
                    //}

                    if (wigWriter != null) {
                        wigWriter.addData(chr, bucketStartPosition, bucketEndPosition, buffer[0]);
                    }

                    if (coverageHistogram != null) {
                        int[] baseCounts = counter.getBaseCount();
                        for (int i = 0; i < baseCounts.length; i++) {
                            coverageHistogram.addDataPoint(baseCounts[i]);
                        }
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
        float nomateCount = 0;
        float posOrientationCount = 0;
        float negOrientationCount = 0;
        String chr;
        int start;
        int end;
        byte[] ref;
        int[] indelCount;
        int[] baseCount;
        int[] baseMismatchCount;

        int ppCount = 0;
        int isizeCount = 0;
        float isizeZScore = 0;
        float totalIsize = 0;

        Counter(String chr, int start, int end) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            baseCount = new int[end - start];
            indelCount = new int[end - start];
            baseMismatchCount = new int[end - start];

            // TODO -- obviously remove this hardcoded genomem reference
            ref = SequenceManager.readSequence("hg18", chr, start, end);
        }

        void increment(int position, byte base) {
            int offset = position - start;
            byte refBase = ref[offset];
            getBaseCount()[offset]++;
            if (refBase != base) {
                baseMismatchCount[offset]++;
            }
            count++;
        }

        void incrementPairedCount() {
            ppCount++;
        }

        void incrementIndel(int position) {
            int offset = position - start;
            indelCount[offset]++;
        }

        void incrementOrientation(int orientation) {
            if (orientation > 0) {
                posOrientationCount++;
            } else if (orientation < 0) {
                negOrientationCount++;
            }
        }


        public void incrementNomate() {
            nomateCount++;
        }

        float getNomateFraction() {
            return (ppCount == 0 ? 0 : nomateCount / ppCount);
        }

        int getNegCount() {
            return negCount;
        }

        public void incrementISize(int inferredInsertSize) {
            //float zs = Math.min(6, (Math.abs(inferredInsertSize) - meanInsertSize) / stdDevInsertSize);
            //isizeZScore += zs;
            //ppCount++;
            //if (Math.abs(Math.abs(inferredInsertSize) - meanInsertSize) > 2 * stdDevInsertSize) {
            //    isizeCount++;
            //}
            totalIsize += Math.abs(inferredInsertSize);
        }

        void incrementNeg() {
            negCount++;
        }

        int getCount() {
            return count;
        }


        float getMaxMismatchFraction() {
            float max = 0.0f;
            for (int i = 0; i < baseMismatchCount.length; i++) {
                max = Math.max(max, (float) baseMismatchCount[i] / getBaseCount()[i]);
            }
            return max;
        }

        float getMaxIndelFraction() {
            float max = 0.0f;
            for (int i = 0; i < indelCount.length; i++) {
                max = Math.max(max, (float) indelCount[i] / getBaseCount()[i]);
            }
            return max;
        }

        float getISizeFraction() {
            if (ppCount < 3) {
                return 0;
            }
            float frac = ((float) isizeCount) / ppCount;
            return frac;
            //float avg = isizeZScore / ppCount;
            //return avg;
        }

        float getAvgIsize() {
            return ppCount == 0 ? 0 : totalIsize / ppCount;
        }

        public int[] getBaseCount() {
            return baseCount;
        }

        public float getPosOrientationFraction() {
            return ppCount == 0 ? 0 : posOrientationCount / ppCount;
        }

        public float getNegOrientationFraction() {
            return ppCount == 0 ? 0 : negOrientationCount / ppCount;
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
        boolean keepZeroes = false;

        WigWriter(File file, int step, boolean keepZeroes) throws IOException {
            this.keepZeroes = keepZeroes;
            this.step = step;
            this.span = step;
            pw = new PrintWriter(new FileWriter(file));
        }

        public void addData(String chr, int start, int end, float data) {

            if (Float.isNaN(data)) {
                return;
            }
            if (genome.getChromosome(chr) == null) {
                return;
            }

            if ((!keepZeroes && data == 0) || end <= start) {
                return;
            }

            int dataSpan = end - start;

            if (chr == null || !chr.equals(lastChr) || dataSpan != span) {
                span = dataSpan;
                outputStepLine(chr, start + 1);
            }
            pw.println((start + 1) + "\t" + data);
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


    // Temporary hack,  illumina only 1 = duplication/translocation,  0 = normal, -1 = inverstion

    private int getOrientation(String oStr) {

        if (oStr.equals("R1F2") || oStr.equals("R2F1")) {
            return 1;
        }
        if (oStr.equals("F1F2") || oStr.equals("F2F1") ||
                oStr.equals("R1R2") || oStr.equals("R2R1")) {
            return -1;
        }
        return 0;
    }


}
