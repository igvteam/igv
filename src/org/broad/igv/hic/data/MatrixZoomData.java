package org.broad.igv.hic.data;

import org.broad.igv.hic.tools.HiCTools;
import org.broad.igv.hic.tools.Preprocessor;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 * @date Aug 10, 2010
 */
public class MatrixZoomData {

    private int chr1;  // Redundant, but convenient    BinDatasetReader
    private int chr2;  // Redundant, but convenient

    private int zoom;
    private int binSize;         // bin size in bp
    private int blockBinCount;   // block size in bins
    private int blockColumnCount;     // number of block columns

    // TODO -- isnt this a memory leak?  Should these be stored?
    private LinkedHashMap<Integer, Block> blocks;


    private Map<Integer, Preprocessor.IndexEntry> blockIndex;
    private DatasetReader reader;


    /**
     * Constructor used by the binary data reader.
     *
     * @param chr1
     * @param chr2
     * @param binSize  in bp
     * @param blockBinCount  block size in bins
     * @param blockColumnCount
     * @param zoom
     * @param blockIndex
     * @param reader
     */
    public MatrixZoomData(int chr1, int chr2, int binSize, int blockBinCount, int blockColumnCount, int zoom,
                          Map<Integer, Preprocessor.IndexEntry> blockIndex, DatasetReader reader) {

        this.chr1 = chr1;
        this.chr2 = chr2;
        this.binSize = binSize;
        this.blockColumnCount = blockColumnCount;
        this.zoom = zoom;
        this.blockBinCount = blockBinCount;
        this.blockIndex = blockIndex;
        blocks = new LinkedHashMap<Integer,Block>(blockIndex.size());
        this.reader = reader;
    }

    public int getBinSize() {
        return binSize;
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

    public int getBlockBinCount() {
        return blockBinCount;
    }

    public int getBlockColumnCount() {
        return blockColumnCount;
    }

    public Map<Integer, Block> getBlocks() {
        return blocks;
    }

    /**
     * Return the blocks overlapping the rectangular region specified.  The units are "bins"
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public List<Block> getBlocksOverlapping(int x1, int y1, int x2, int y2) {

        int col1 = x1 / blockBinCount;
        int row1 = y1 / blockBinCount;

        int col2 = x2 / blockBinCount;
        int row2 = y2 / blockBinCount;

        int maxSize = (col2 - col1 + 1) * (row2 - row1 + 1);

        List<Block> blockList = new ArrayList<Block>(maxSize);
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
            if (reader != null && blockIndex != null) {
                b = readBlock(blockNumber);
                blocks.put(blockNumber, b);
            }
        }
        return b;
    }

    private Block readBlock(int blockNumber) {
        Preprocessor.IndexEntry idx = blockIndex.get(blockNumber);
        Block b = null;
        if (idx != null) {
            try {
                b = reader.readBlock(blockNumber, idx);

            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        return b;
    }


    // Dump the contents to standard out
    public void dump(Chromosome[] chromosomes) {

        // Get the block index keys, and sort
        List<Integer> blockNumbers = new ArrayList<Integer>(blockIndex.keySet());
        Collections.sort(blockNumbers);

        System.out.println("# " + chromosomes[chr1].getName() + " - " + chromosomes[chr2].getName());

        for(int blockNumber : blockNumbers) {
            Block b = readBlock(blockNumber);
            if(b != null) {
                for(ContactRecord rec : b.getContactRecords()) {
                     System.out.println(rec.getX()*binSize + "\t" + rec.getY()*binSize + "\t" + rec.getCounts());
                }
            }
        }

    }
}
