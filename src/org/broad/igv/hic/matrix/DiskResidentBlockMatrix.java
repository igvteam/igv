package org.broad.igv.hic.matrix;

import org.broad.igv.hic.MainWindow;
import org.broad.igv.util.ObjectCache;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.io.*;

/**
 * Field	Type	Description
 * Magic number	Integer	Value should be 6515048, which when read as a zero delimited string is “hic”.
 * Version	Integer	The version number, currently =  1
 * Genome	String	The genome ID (e.g. “hg19”)
 * Chromosome 1	String	Name of first chromosome
 * Chromosome 2	String	Name of second chromosome
 * Bin size	Integer	The bin size in base-pairs
 * Lower bounds for scale	Float	5th percentile suggested
 * Upper bounds for scale	Float	95th percentile suggested
 * Number of rows	Integer	Number of rows in the matrix.  Rows correspond to chromosome 1.
 * Number of columns	Integer
 * Block size  Integer
 *
 * @author jrobinso
 *         Date: 8/8/12
 *         Time: 8:42 AM
 */
public class DiskResidentBlockMatrix implements BasicMatrix {

    String path;
    private String genome;
    private String chr1;
    private String chr2;
    private int binSize;
    float lowerValue;
    float upperValue;
    int dim;
    int blockSize;
    int remSize;   // Dimension of last block

    int arrayStartPosition;
    boolean isLoading = false;

    ObjectCache<String, float[][]> blockDataCache = new ObjectCache<String, float[][]>(200);
    private int nFullBlocks;

    public DiskResidentBlockMatrix(String path) throws IOException {
        this.path = path;
        init();
    }

    public String getChr1() {
        return chr1;
    }

    @Override
    public float getEntry(int row, int col) {

        int blockRowIdx = row / blockSize;
        int blockColIdx = col / blockSize;
        String key = "row" + blockRowIdx + "_col" + blockColIdx;
        float[][] blockData = blockDataCache.get(key);
        if (blockData == null) {
            MainWindow.getInstance().showGlassPane();
            blockData = loadBlockData(blockRowIdx, blockColIdx);
            blockDataCache.put(key, blockData);
            MainWindow.getInstance().hideGlassPane();
        }

        if (blockData == null) {
            return Float.NaN;
        } else {
            int rowRelative = row - blockRowIdx * blockSize;
            int colRelative = col - blockColIdx * blockSize;
            return blockData[rowRelative][colRelative];
        }

    }

    public synchronized float[][] loadBlockData(int blockRowIdx, int blockColIdx) {

        String key = "row" + blockRowIdx + "_col" + blockColIdx;
        float [][] blockData = blockDataCache.get(key);
        if(blockData != null) return blockData;    // In case this was calculated in another thread

        SeekableStream is = null;
        try {
            is = SeekableStreamFactory.getStreamFor(path);

            int pointsPerBlockRow = blockSize * dim;  // Applies to all but the last row

            int rowDim = blockRowIdx < nFullBlocks ? blockSize : remSize;
            int colDim = blockColIdx < nFullBlocks ? blockSize : remSize;

            long startFilePosition =
                    arrayStartPosition + (blockRowIdx * pointsPerBlockRow + blockColIdx * blockSize * rowDim) * 4;


             int nDataPoints = rowDim * colDim;
            int nBytes = nDataPoints * 4;
            byte[] byteArray = new byte[nBytes];

            System.out.println("Loading block " + blockRowIdx + " " + blockColIdx + "  rowDim=" + rowDim +
                    "  colDim=" + colDim + " nDataPoints=" + nDataPoints + "  startFilePosition=" + startFilePosition
            + "  endFilePosition=" + (startFilePosition + nBytes) +
            "  thread=" + Thread.currentThread().getName());

            is.seek(startFilePosition);
            is.readFully(byteArray);

            ByteArrayInputStream bis = new ByteArrayInputStream(byteArray);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);


            blockData = new float[rowDim][colDim];

            for (int r = 0; r < rowDim; r++) {
                for (int c = 0; c < colDim; c++) {
                    float f = les.readFloat();
                    blockData[r][c] = f;
                }
            }

            blockDataCache.put(key, blockData);
            return blockData;

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            return null;
        } finally {
            if (is != null)
                try {
                    is.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
        }
    }


    @Override
    public int getRowDimension() {
        return dim;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getColumnDimension() {
        return dim;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public float getLowerValue() {
        return lowerValue;
    }

    public float getUpperValue() {
        return upperValue;
    }


    @Override
    public BasicMatrix getSubMatrix(int startRow, int endRow, int startCol, int endCol) {
        return null;
    }

    public int getRemSize() {
        return remSize;
    }

    private void init() throws IOException {
        // Peak at file to determine version
        BufferedInputStream bis = null;
        try {
            InputStream is = ParsingUtils.openInputStream(path);
            bis = new BufferedInputStream(is);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);

            int bytePosition = 0;
            int magic = les.readInt();    // <= 6515048
            bytePosition += 4;
            System.out.println("Magic number = " + magic);

            //if (magic != 6515048)      <= ERROR
            // Version number
            int version = les.readInt();
            bytePosition += 4;
            System.out.println("Version = " + version);

            genome = les.readString();
            bytePosition += genome.length() + 1;
            System.out.println("Genome = " + genome);

            chr1 = les.readString();
            bytePosition += chr1.length() + 1;
            System.out.println("Chr1 = " + chr1);

            chr2 = les.readString();
            bytePosition += chr2.length() + 1;
            System.out.println("Chr2 = " + chr2);

            binSize = les.readInt();
            bytePosition += 4;
            System.out.println("binSize = " + binSize);

            lowerValue = les.readFloat();
            bytePosition += 4;
            System.out.println("Lower value = " + lowerValue);

            upperValue = les.readFloat();
            bytePosition += 4;
            System.out.println("Upper value = " + upperValue);

            int nRows = les.readInt();  // # rows, assuming square matrix
            bytePosition += 4;
            System.out.println("Row count = " + nRows);

            int nCols = les.readInt();
            bytePosition += 4;
            System.out.println("Column count = " + nRows);

            if (nRows == nCols) {
                dim = nRows;
            } else {
                throw new RuntimeException("Non-square matrices not supported");
            }

            blockSize = les.readInt();
            bytePosition += 4;
            System.out.println("Block size = " + blockSize);

            nFullBlocks = dim / blockSize;
            remSize = dim - nFullBlocks * blockSize;
            System.out.println(nFullBlocks + "  " + remSize);

            this.arrayStartPosition = bytePosition;


        } finally {
            if (bis != null)
                bis.close();
        }
    }


    public static void main(String [] args) throws IOException {
        DiskResidentRowMatrix bm = new DiskResidentRowMatrix("/Users/jrobinso/projects/hic/pearsons_14__14_250000.bin");
        saveAsBlockMatrix(bm, 10, new File("/Users/jrobinso/projects/hic/block_test_chr14_1M.bin"));
    }

    public static void saveAsBlockMatrix(DiskResidentRowMatrix bm, int blockSize, File outputFile) throws IOException {

        LittleEndianOutputStream los = null;
        try {
            BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(outputFile));
            los = new LittleEndianOutputStream(bos);

            los.writeInt(6515048);    // magic number
            los.writeInt(2);          // version


            los.writeString(bm.getGenome());
            los.writeString(bm.getChr1());
            los.writeString(bm.getChr2());
            los.writeInt(bm.getBinSize());
            los.writeFloat(bm.getLowerValue());
            los.writeFloat(bm.getUpperValue());
            los.writeInt(bm.getDim());
            los.writeInt(bm.getDim());
            los.writeInt(blockSize);

            int nFullBlocks = bm.getDim() / blockSize;
            int remSize = bm.getDim() - nFullBlocks * blockSize;

            // Loop through blocks
            for (int blockRowIdx = 0; blockRowIdx <= nFullBlocks; blockRowIdx++) {

                int startRow = blockRowIdx * blockSize;
                int rowDim = blockRowIdx < nFullBlocks ? blockSize : remSize;

                for (int blockColIdx = 0; blockColIdx <= nFullBlocks; blockColIdx++) {

                    int startCol = blockColIdx * blockSize;
                    int colDim = blockColIdx < nFullBlocks ? blockSize : remSize;

                    for (int r = 0; r < rowDim; r++) {
                        for (int c = 0; c < colDim; c++) {
                            int row = startRow + r;
                            int col = startCol + c;

                            float value = bm.getEntry(row, col);
                            los.writeFloat(value);
                        }
                    }
                }
            }


        } finally {
            if (los != null)
                los.close();
        }
    }

}
