package org.broad.igv.hic.data;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.hic.HiC;
import org.broad.igv.hic.matrix.BasicMatrix;
import org.broad.igv.hic.matrix.DiskResidentBlockMatrix;
import org.broad.igv.hic.matrix.DiskResidentRowMatrix;
import org.broad.igv.hic.matrix.InMemoryMatrix;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.collections.DoubleArrayList;
import org.broad.igv.util.collections.DownsampledDoubleArrayList;
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
        //String path = "/Users/jrobinso/projects/hic/bin_chr14_1M.bin";
        String path = "/Users/jrobinso/projects/hic/chr14_5e3_N17658_output 2.bin";
        //String path = "/Users/jrobinso/projects/hic/pearsons_14__14_100000.bin";
        //createBlockIndexedFile(f, null, 50);
        //BasicMatrix bm = readPearsons(path);
        // System.out.println(bm.getColumnDimension());
        //readPearsons(path);
        //writeHeader("/Users/jrobinso/projects/hic/header.bin");
        testReadPy();
    }


    public static void testReadPy() throws IOException {
        // Peak at file to determine version
        BufferedInputStream bis = null;

        InputStream is = ParsingUtils.openInputStream("/Users/jrobinso/foo.txt");
        bis = new BufferedInputStream(is);
        LittleEndianInputStream les = new LittleEndianInputStream(bis);

        int b;
        while ((b = les.read()) >= 0) {
            System.out.println(b + "  " + ((char) b));
        }

//        String s = les.readString();
//        System.out.println(s);
        is.close();
    }

    public static BasicMatrix readPearsons(String path) throws IOException {

        // Peak at file to determine version
        BufferedInputStream bis = null;
        int magic;
        int version;
        try {
            InputStream is = ParsingUtils.openInputStream(path);
            bis = new BufferedInputStream(is);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);

            magic = les.readInt();

            if (magic == 6515048) {
                version = les.readInt();
            } else {
                throw new RuntimeException("Unsupported format: " + path);
            }
        } finally {
            if (bis != null)
                bis.close();
        }

        if (version == 1) {
            return new DiskResidentRowMatrix(path);

        } else {
            return new DiskResidentBlockMatrix(path);
        }


    }


    public static void estimatePercentiles(String path) throws IOException {
        BufferedInputStream bis = null;
        float max = -1;
        float min = 1;
        DownsampledDoubleArrayList dataSampleList = new DownsampledDoubleArrayList(100000);

        try {
            InputStream is = ParsingUtils.openInputStream(path);
            bis = new BufferedInputStream(is);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);

            int dim = les.readInt();
            int total = dim * dim;
            for (int i = 0; i < total; i++) {
                float f = les.readFloat();
                if (f < 1) {
                    max = f > max ? f : max;
                    min = f < min ? f : min;
                    dataSampleList.add((double) f);
                }
            }

            double[] data = dataSampleList.toArray();
            System.out.println("min = " + min);
            System.out.println("2.5 % = " + StatUtils.percentile(data, 2.5));
            System.out.println("5 % = " + StatUtils.percentile(data, 5));
            System.out.println("10 % = " + StatUtils.percentile(data, 10));
            System.out.println("90 % = " + StatUtils.percentile(data, 90));
            System.out.println("95 % = " + StatUtils.percentile(data, 95));
            System.out.println("97.5 % = " + StatUtils.percentile(data, 97.5));
            System.out.println("max = " + max);

            bis.close();

        } finally {
            if (bis != null) bis.close();

        }
    }


    public static void convertYunfanFormat(String path, String chr, int binSize) throws IOException {
        BufferedInputStream bis = null;

        try {
            InputStream is = ParsingUtils.openInputStream(path);
            bis = new BufferedInputStream(is);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);

            int dim = les.readInt();
            int nPoints = dim * dim;
            float[] data = new float[nPoints];
            for (int i = 0; i < nPoints; i++) {
                data[i] = les.readFloat();
            }

            BasicMatrix bm = new InMemoryMatrix(data, dim);
            File f = new File("pearsons_" + chr + "_" + chr + "_" + binSize);
            dumpPearsonsBinary(bm, chr, chr, binSize, f);

        } finally {
            if (bis != null) bis.close();

        }
    }


    public static void writeHeader(String path) throws IOException {

        File f = new File(path);

        FileOutputStream fos = new FileOutputStream(f);
        BufferedOutputStream bos = new BufferedOutputStream(fos);
        LittleEndianOutputStream los = new LittleEndianOutputStream(bos);

        // Magic number - 4 bytes
        los.writeByte('h');
        los.writeByte('i');
        los.writeByte('c');
        los.writeByte(0);

        // Version number
        los.writeInt(1);

        // Genome --
        los.writeString("hg19");

        // Chromosomes
        los.writeString("14");
        los.writeString("14");

        // Resolution (bin size)
        los.writeInt(5000);

        // Statistics, other attributes
        los.writeFloat(-0.004103539418429137f);
        los.writeFloat(0.03536746241152287f);
        los.writeInt(21458);  // # rows, assuming square matrix

        los.close();
        bos.close();
        fos.close();

    }

    /**
     * Dump the pearson's correlation -- for development
     *
     * @param pearsons
     * @param chr1
     * @param chr2
     */
    public static void dumpPearsonsBinary(BasicMatrix pearsons, String chr1, String chr2, int binSize, File f) throws IOException {

        int nCols = pearsons.getColumnDimension();
        // Assuming sqaure matrix for this test

        OutputStream fos = null;

        fos = new FileOutputStream(f);
        BufferedOutputStream bos = new BufferedOutputStream(fos);
        LittleEndianOutputStream los = new LittleEndianOutputStream(bos);

        // Magic number - 4 bytes
        los.writeByte('h');
        los.writeByte('i');
        los.writeByte('c');
        los.writeByte(0);

        // Version number
        los.writeInt(1);

        // Genome --
        los.writeString("hg19");

        // Chromosomes
        los.writeString(chr1);
        los.writeString(chr2);

        // Resolution (bin size)
        los.writeInt(binSize);

        // Statistics, other attributes
        DoubleArrayList flattenedDataList = new DoubleArrayList(pearsons.getColumnDimension() * pearsons.getRowDimension());
        double min = 1;
        double max = -1;
        for (int i = 0; i < pearsons.getRowDimension(); i++) {
            for (int j = 0; j < pearsons.getColumnDimension(); j++) {
                double value = pearsons.getEntry(i, j);
                if (!Double.isNaN(value) && value != 1) {
                    min = value < min ? value : min;
                    max = value > max ? value : max;
                    flattenedDataList.add(value);
                }
            }
        }

        // Stats
        double[] flattenedData = flattenedDataList.toArray();
        los.writeFloat((float) StatUtils.percentile(flattenedData, 5));
        los.writeFloat((float) StatUtils.percentile(flattenedData, 95));


        // Data
        los.writeInt(nCols);  // # rows, assuming square matrix
        los.writeInt(nCols);
        for (int r = 0; r < pearsons.getRowDimension(); r++) {
            for (int c = 0; c < pearsons.getColumnDimension(); c++) {
                los.writeFloat(pearsons.getEntry(r, c));
            }
        }
        los.close();
        bos.close();
        fos.close();

    }


    /**
     * Dump the eigenvector  correlation -- for development
     */
    public static void dumpEigenvector(HiC hic) throws IOException {

        int step = hic.zd.getBinSize();
        int chrIdx = hic.zd.getChr1();
        String chr = hic.getChromosomes()[chrIdx].getName();
        double[] eigenvector = hic.getEigenvector(0);
        // Assuming sqaure matrix for this test

        OutputStream fos = null;

        String fn = "eigenvector_" + chr + "_" + step + ".wig";


        PrintWriter pw = new PrintWriter(new FileWriter(fn));
        pw.println("track type=wiggle name=Eigenvector");

        int start = 0;
        // Skip throuh NaN
        int idx = 0;
        for (; idx < eigenvector.length; idx++) {
            if (Double.isNaN(eigenvector[idx])) {
                start += step;
            } else {
                break;
            }
        }

        pw.println("fixedStep chrom=" + chr + " start=" + start + " step=" + step + " span=" + step);

        for (; idx < eigenvector.length; idx++) {
            pw.println(eigenvector[idx]);
        }

        pw.close();

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
