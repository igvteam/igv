package org.igv.sam;

import org.apache.commons.math3.stat.StatUtils;
import org.igv.util.collections.DoubleArrayList;

import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Created by jrobinso on 12/20/16.
 * <p>
 * Some simple statistics to distinguish dna, rna, and 3rd gen alignment datasets.
 */


public class ReadStats {

    private static int MAX_READ_COUNT = 100000;   // Sufficiently high to get good statistics

    public int readCount = 0;
    private int indelCount = 0;
    private int nCount = 0;
    private DoubleArrayList readLengths = new DoubleArrayList(10000);
    private DoubleArrayList refToReadRatios = new DoubleArrayList(10000);
    private Set<Integer> qualityScores = new java.util.HashSet();

    public double medianReadLength = 0;
    public double readLengthStdDev = 0;
    public double medianRefToReadRatio = 0;
    public double fracReadsWithIndels;
    public double fracReadsWithNs;
    public double averageCigarLength;

    private static final Random RAND = new Random();

    public void addAlignment(Alignment alignment) {

        if (readCount > MAX_READ_COUNT) return;

        if (alignment.isMapped()) {

            readCount++;

            final int readLength = alignment.getReadSequence().length();
            final char[] cigarString = alignment.getCigarString().toCharArray();

            readLengths.add(readLength);

            if (containsChar(cigarString, 'N')) {
                nCount++;
                refToReadRatios.add(((double) (alignment.getAlignmentEnd() - alignment.getAlignmentStart())) / readLength);
            }

            if (containsChar(cigarString, 'D') || containsChar(cigarString, 'I')) {
                indelCount++;
            }

            double cigarLength = alignment.getCigarString().length();
            averageCigarLength = averageCigarLength * (((double) readCount - 1) / readCount) + (cigarLength / readCount);
        }

        for(AlignmentBlock ab : alignment.getAlignmentBlocks()) {
            if(ab.hasQualities()) {
                for (int i = 0; i < ab.getLength(); i++) {
                    qualityScores.add((int) ab.getQuality(i));
                }
            }
        }
    }

    public void compute() {

        if (readLengths.size() > 0) {
            final double[] downsampled = downsample(readLengths, 10000);
            medianReadLength = StatUtils.percentile(downsampled, 50);
            readLengthStdDev = Math.sqrt(StatUtils.variance(downsampled));
        }


        if (refToReadRatios.size() > 0) {
            medianRefToReadRatio = StatUtils.percentile(downsample(refToReadRatios, 10000), 50);
        }

        fracReadsWithIndels = ((double) indelCount) / readCount;
        fracReadsWithNs = ((double) nCount) / readCount;

    }

    // Hand coded for speed
    private boolean containsChar(char[] charArray, char c) {
        for (int i = 0; i < charArray.length; i++) {
            if (charArray[i] == c) return true;
        }
        return false;
    }

    private double[] downsample(DoubleArrayList list, int size) {

        if (list.size() < size) return list.toArray();

        else {
            double[] ds = list.toArray(0, size);

            for (int i = size; i < list.size(); i++) {
                double samplingProb = ((double) size) / (size + (i - size) + 1);
                if (RAND.nextDouble() < samplingProb) {
                    int idx = (int) (RAND.nextDouble() * (size - 1));
                    ds[idx] = list.get(i);

                }
            }

            return ds;
        }

    }

    public AlignmentTrack.ExperimentType inferType() {

        // SBX determination does not require statistics.  Try this first.
        if(qualityScores.size() > 0 && SbxUtils.isSBXAlignments(qualityScores)) {
            return AlignmentTrack.ExperimentType.SBX;
        }

        compute();
        if (readCount < 20) return AlignmentTrack.ExperimentType.UNKOWN; // Not enough reads
        if (medianRefToReadRatio > 5 || fracReadsWithNs > 0.01) {
            return AlignmentTrack.ExperimentType.RNA;
        } else if ((readLengthStdDev > 100 || medianReadLength > 1000) && averageCigarLength > 10) { // Cigar length to filter consensus reads
            return AlignmentTrack.ExperimentType.THIRD_GEN;
        } else {
            return AlignmentTrack.ExperimentType.OTHER;
        }

    }

}
