package org.igv.hic;

import org.igv.feature.Chromosome;
import org.igv.ucsc.twobit.UnsignedByteBuffer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents zoom-level data for a matrix.
 * Assumes Chromosome has methods: getName(), getSize().
 * Assumes Region has methods: getChr(), getStart(), getEnd().
 * Assumes UnsignedByteBuffer has methods: getString(), getInt(), getFloat(), getLong().
 */
public class MatrixZoomData {

    private final Chromosome chr1;
    private final Chromosome chr2;

    // zoom info
    public static record Zoom(int index, String unit, int binSize) {
    }

    private Zoom zoom;
    int blockBinCount;
    int blockColumnCount;
    private StaticBlockIndex blockIndex;

    double averageCount;
    double sumCounts;
    double stdDev;
    double occupiedCellCount;
    double percent95;

    public MatrixZoomData(Chromosome chr1, Chromosome chr2) {
        this.chr1 = chr1;
        this.chr2 = chr2;
    }

    public String getKey() {
        return chr1.getName() + "_" + chr2.getName() + "_" + zoom.unit + "_" + zoom.binSize;
    }

    public List<Integer> getBlockNumbers(Region region1, Region region2, int version) {
        // swap if regions are reversed
        if (region1.chr().equals(this.chr2.getName()) && region2.chr().equals(this.chr1.getName())) {
            Region tmp = region1;
            region1 = region2;
            region2 = tmp;
        }
        boolean sameChr = this.chr1 == this.chr2;
        int binsize = zoom.binSize;
        int blockBinCount = this.blockBinCount;
        int blockColumnCount = this.blockColumnCount;

        if (version < 9 || !sameChr) {
            return getBlockNumbersV8(region1, region2, binsize, blockBinCount, blockColumnCount, sameChr);
        } else {
            return getBlockNumbersV9(region1, region2, binsize, blockBinCount, blockColumnCount);
        }
    }

    private List<Integer> getBlockNumbersV8(Region region1, Region region2,
                                            int binsize, int blockBinCount, int blockColumnCount, boolean sameChr) {
        long x1 = region1.start() / binsize;
        long x2 = region1.end() / binsize;
        long y1 = region2.start() / binsize;
        long y2 = region2.end() / binsize;

        int col1 = (int) Math.floor((double) x1 / blockBinCount);
        int col2 = (int) Math.floor((double) (x2 - 1) / blockBinCount);
        int row1 = (int) Math.floor((double) y1 / blockBinCount);
        int row2 = (int) Math.floor((double) (y2 - 1) / blockBinCount);

        List<Integer> blockNumbers = new ArrayList<>();
        for (int row = row1; row <= row2; row++) {
            for (int column = col1; column <= col2; column++) {
                int blockNumber;
                if (sameChr && row < column) {
                    blockNumber = column * blockColumnCount + row;
                } else {
                    blockNumber = row * blockColumnCount + column;
                }
                if (!blockNumbers.contains(blockNumber)) {
                    blockNumbers.add(blockNumber);
                }
            }
        }
        return blockNumbers;
    }

    private List<Integer> getBlockNumbersV9(Region region1, Region region2,
                                            int binsize, int blockBinCount, int blockColumnCount) {
        long binX1 = region1.start() / binsize;
        long binX2 = region1.end() / binsize;
        long binY1 = region2.start() / binsize;
        long binY2 = region2.end() / binsize;

        int translatedLowerPAD = (int) Math.floor((double) (binX1 + binY1) / 2 / blockBinCount);
        int translatedHigherPAD = (int) Math.floor((double) (binX2 + binY2) / 2 / blockBinCount);
        int translatedNearerDepth = (int) Math.floor(Math.log(1 + Math.abs(binX1 - binY2) / Math.sqrt(2) / blockBinCount) / Math.log(2));
        int translatedFurtherDepth = (int) Math.floor(Math.log(1 + Math.abs(binX2 - binY1) / Math.sqrt(2) / blockBinCount) / Math.log(2));

        boolean containsDiagonal = (binX2 - binY1) * (binX1 - binY2) < 0;
        int nearerDepth = containsDiagonal ? 0 : Math.min(translatedNearerDepth, translatedFurtherDepth);
        int furtherDepth = Math.max(translatedNearerDepth, translatedFurtherDepth);

        List<Integer> blockNumbers = new ArrayList<>();
        for (int depth = nearerDepth; depth <= furtherDepth; depth++) {
            for (int pad = translatedLowerPAD; pad <= translatedHigherPAD; pad++) {
                int blockNumber = depth * blockColumnCount + pad;
                blockNumbers.add(blockNumber);
            }
        }
        return blockNumbers;
    }

    public static MatrixZoomData parseMatrixZoomData(Chromosome chr1, Chromosome chr2, UnsignedByteBuffer dis) throws IOException {
        MatrixZoomData zd = new MatrixZoomData(chr1, chr2);

        String unit = dis.getString();
        int zoomIndex = dis.getInt();
        double sumCounts = dis.getFloat();
        double occupiedCellCount = dis.getFloat();
        double stdDev = dis.getFloat();
        double percent95 = dis.getFloat();
        int binSize = dis.getInt();
        zd.blockBinCount = dis.getInt();
        zd.blockColumnCount = dis.getInt();
        int nBlocks = dis.getInt();

        zd.zoom = new Zoom(zoomIndex, unit, binSize);
        zd.blockIndex = new StaticBlockIndex(nBlocks, dis);

        int nBins1 = (chr1.getLength() / binSize);
        int nBins2 = (chr2.getLength() / binSize);
        double avgCount = (sumCounts / nBins1) / nBins2;

        zd.averageCount = avgCount;
        zd.sumCounts = sumCounts;
        zd.stdDev = stdDev;
        zd.occupiedCellCount = occupiedCellCount;
        zd.percent95 = percent95;

        return zd;
    }

    // getters if needed
    public Chromosome getChr1() {
        return chr1;
    }

    public Chromosome getChr2() {
        return chr2;
    }

    public Zoom getZoom() {
        return zoom;
    }

    public StaticBlockIndex getBlockIndex() {
        return blockIndex;
    }

    public double getAverageCount() {
        return averageCount;
    }

    public double getSumCounts() {
        return sumCounts;
    }
}