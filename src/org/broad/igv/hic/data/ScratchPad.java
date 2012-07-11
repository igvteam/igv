package org.broad.igv.hic.data;

import org.apache.commons.math.linear.RealMatrix;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;

/**
 * @author jrobinso
 *         Date: 6/29/12
 *         Time: 2:27 PM
 */
public class ScratchPad {

    public static void main(String[] args) throws IOException {
        File f = new File("/Users/jrobinso/projects/hic/bin_chr14_1M.bin");
        //File f = new File("/Users/jrobinso/projects/hic/chr14_5e3_N17658_output.bin");
        //File f = new File("/Users/jrobinso/projects/hic/chr14_50K_test.bin");
        createBlockIndexedFile(f, null, 50);
    }

    public static void readPearsons(File file) throws IOException {

        FileInputStream fis = null;

        fis = new FileInputStream(file);
        BufferedInputStream bis = new BufferedInputStream(fis);
        LittleEndianInputStream les = new LittleEndianInputStream(bis);

        int nRows = les.readInt();
        int nTot = nRows * nRows;
        int nBytes = 4;
        for (int i = 0; i < nTot; i++) {
            //les.readByte();
            Float f = les.readFloat();
            System.out.println(f);
            nBytes += 4;

        }
        System.out.println("Nbytes=" + nBytes);


        fis.close();
    }

    /**
     * Dump the pearson's correlation -- for development
     *
     * @param pearsons
     */
    public static void dumpPearsonsBinary(RealMatrix pearsons) throws IOException {

        int nCols = pearsons.getColumnDimension();
        // Assuming sqaure matrix for this test

        OutputStream fos = null;

        fos = new FileOutputStream("test_out.bin");
        BufferedOutputStream bos = new BufferedOutputStream(fos);
        LittleEndianOutputStream los = new LittleEndianOutputStream(bos);

        los.writeInt(nCols);

        double[][] data = pearsons.getData();
        for (int r = 0; r < pearsons.getRowDimension(); r++) {
            for (int c = 0; c < pearsons.getColumnDimension(); c++) {
                los.writeFloat((float) data[r][c]);
            }
        }
        los.close();
        bos.close();
        fos.close();

    }


    public static void createBlockIndexedFile(File inputFile, File outputFile, int blockSize) throws IOException {

        FileInputStream fis = null;

        fis = new FileInputStream(inputFile);
        BufferedInputStream bis = new BufferedInputStream(fis);
        LittleEndianInputStream les = new LittleEndianInputStream(bis);

        int nDataRows = les.readInt();
        int nDataColumns = nDataRows;
        int nTot = nDataRows * nDataColumns;

        int nBlockColumns = nDataColumns / blockSize + 1;
        int nBlockRows = nDataRows / blockSize + 1;

        // The end blocks will in general be smaller
        int lastRowBlockSize = nDataRows - (nBlockRows - 1) * blockSize;
        int lastColBlockSize = nDataColumns - (nBlockColumns - 1) * blockSize;

        // Now transform matrix from linear -> block layout
        // Data is read in row-major order,  tranforming rows of blocks simultaneously

        int lastBlockRow = 0;

        int nBlockRow = 0;
        Block[] blocksForRow = new Block[nBlockColumns];
        resetBlocks(blocksForRow, blockSize, nBlockColumns, lastColBlockSize, nBlockRow);

        for (int i = 0; i < nTot; i++) {
            //les.readByte();

            int row = i / nDataColumns;
            nBlockRow = row / blockSize;

            int col = (i - row * nDataColumns);
            int nBlockCol = col / blockSize;

            if (nBlockRow != lastBlockRow) {
                for (Block b : blocksForRow) {
                    System.out.println(lastBlockRow + "\t" + nBlockCol + "\t" + b.data.length);
                }
                resetBlocks(blocksForRow, blockSize, nBlockColumns, lastColBlockSize, nBlockRow);
                lastBlockRow = nBlockRow;
            }

            // Index within block
            int startIdx = nBlockRow * blockSize * nDataColumns + nBlockCol * blockSize;
            int idx = i - startIdx;

            Float f = les.readFloat();
            blocksForRow[nBlockCol].add(f);
        }

        for (Block b : blocksForRow) {
            System.out.println(nBlockRow + "\t" + "" + "\t" + b.data.length);
        }
        resetBlocks(blocksForRow, blockSize, nBlockColumns, lastColBlockSize, nBlockRow);
        lastBlockRow = nBlockRow;



        fis.close();
    }

    private static void resetBlocks(Block[] blocksForRow,
                                       int blockSize, int nBlockSize, int lastBlockSize, int nBlockRow) {

        int tmp = nBlockRow > nBlockSize - 1 ? lastBlockSize : blockSize;
        for (int i = 0; i < nBlockSize - 1; i++) {
            blocksForRow[i] = new Block(tmp, blockSize);
        }
        blocksForRow[nBlockSize - 1] = new Block(tmp, lastBlockSize);
    }


    static class Block {
        int nRows;
        int nCols;
        float[] data;
        int nextIdx = 0;

        Block(int nRows, int nCols) {
            this.nRows = nRows;
            this.nCols = nCols;
            data = new float[nRows * nCols];
        }

        void add(float d) {
            try {
                data[nextIdx] = d;
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            nextIdx++;
        }
    }
}
