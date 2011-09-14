package org.broad.igv.hic.data;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 * @date Aug 10, 2010
 */
public class MatrixZoomData {

    private int chr1;  // Redundant, but convenient    BinDatasetReader
    private int chr2; // Redundant, but convenient

    private int zoom;
    private int binSize;     // in bp
    private int blockSize;   // in bins
    private int blockColumnCount;   // number of block columns

    private LinkedHashMap<Integer, Block> blocks;
    private Map<Integer, DatasetWriter.IndexEntry> blockIndex;
    private DatasetReader reader ;


    /**
     * Constructor used by the alignment file parser.
     */
    public MatrixZoomData(int chr1, int chr2, int binSize, int blockColumnCount, int zoom) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        this.binSize = binSize;
        this.blockColumnCount = blockColumnCount;
        this.zoom = zoom;

        int nBinsX = Dataset.chromosomes[chr1].getSize() / binSize + 1;
        blockSize = nBinsX / blockColumnCount + 1;
        blocks = new LinkedHashMap(blockColumnCount * blockColumnCount);
    }

    /**
     * Constructor used by the binary data reader.
     * @param chr1
     * @param chr2
     * @param binSize
     * @param blockSize
     * @param blockColumnCount
     * @param zoom
     * @param blockIndex
     * @param reader
     */
    public MatrixZoomData(int chr1, int chr2, int binSize, int blockSize, int blockColumnCount, int zoom,
                          Map<Integer, DatasetWriter.IndexEntry> blockIndex, DatasetReader reader) {

        this.chr1 = chr1;
        this.chr2 = chr2;
        this.binSize = binSize;
        this.blockColumnCount = blockColumnCount;
        this.zoom = zoom;
        this.blockSize = blockSize;
        this.blockIndex = blockIndex;
        blocks = new LinkedHashMap(blockIndex.size());
        this.reader = reader;
    }

    public void incrementCount(int pos1, int pos2) {
        int xBin = pos1 / getBinSize();
        int yBin = pos2 / getBinSize();

        if (chr1 == chr2) {
            int b1 = Math.min(xBin, yBin);
            int b2 = Math.max(xBin, yBin);
            xBin = b1;
            yBin = b2;
        }

        // compute block number (fist block is zero)
        int blockCol = xBin / getBlockSize();
        int blockRow = yBin / getBlockSize();
        int blockNumber = getBlockColumnCount() * blockRow + blockCol;

        Block block = blocks.get(blockNumber);
        if (block == null) {
            block = new Block(blockNumber);
            blocks.put(blockNumber, block);
        }
        block.incrementCount(xBin, yBin);

    }

    public void parsingComplete() {
        for (Block b : blocks.values()) {
            b.parsingComplete();
        }
    }

    public int getBinSize() {
        return binSize;
    }

    public void setBinSize(int binSize) {
        this.binSize = binSize;
    }

    public int getChr1() {
        return chr1;
    }


    public int getChr2() {
        return chr2;
    }

    public int getZoom() {
        return zoom;
    }

    public void setZoom(int zoom) {
        this.zoom = zoom;
    }

    public int getBlockSize() {
        return blockSize;
    }

    public int getBlockColumnCount() {
        return blockColumnCount;
    }

    public Map<Integer, Block> getBlocks() {
        return blocks;
    }


    public List<Block> getBlocksOverlapping(int x1, int y1, int x2, int y2) {

        int col1 = (x1 / binSize) / getBlockSize();
        int row1 = (y1 / binSize) / getBlockSize();

        int col2 = (x2 / binSize) / getBlockSize();
        int row2 = (y2 / binSize) / getBlockSize();

        int maxSize = (col2 - col1 + 1) * (row2 - row1 + 1);

        List<Block> blockList = new ArrayList(maxSize);
        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                int blockNumber = r * getBlockColumnCount() + c;
                Block b = getBlock(blockNumber);
                if (b != null) {
                    blockList.add(b);
                }
            }
        }
        return blockList;
    }

    public Block getBlock(int blockNumber) {
        Block b = blocks.get(blockNumber);
        if (b == null) {
            if(reader != null && blockIndex != null) {
                DatasetWriter.IndexEntry idx = blockIndex.get(blockNumber);
                if(idx != null) {
                    try {
                        b = reader.readBlock(blockNumber, idx);
                        blocks.put(blockNumber, b);
                    } catch (IOException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }
                }
            }
        }
        return b;
    }

}
