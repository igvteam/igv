package org.broad.igv.hic.data;

import org.apache.commons.math.linear.*;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.broad.igv.hic.tools.Preprocessor;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * @author jrobinso
 * @date Aug 10, 2010
 */
public class MatrixZoomData {

    // A thread executor for computing Pearsons correlation
    // TODO -- move this to some utility class


    private Chromosome chr1;  // Redundant, but convenient
    private Chromosome chr2;  // Redundant, but convenient

    private int zoom;
    private int binSize;         // bin size in bp
    private int blockBinCount;   // block size in bins
    private int blockColumnCount;     // number of block columns

    // TODO -- isnt this a memory leak?  Should these be stored?
    private LinkedHashMap<Integer, Block> blocks;


    private Map<Integer, Preprocessor.IndexEntry> blockIndex;
    private DatasetReader reader;
    private RealMatrix pearsons;
    private double pearsonsMin = -1;
    private double pearsonsMax = 1;
    private RealMatrix oe;
    private double[] eigenvector;

    public class ScaleParameters {
        double percentile90;
        double mean;

        ScaleParameters(double mean, double percentile90) {
            this.mean = mean;
            this.percentile90 = percentile90;
        }
    }


    /**
     * Construct from a binary stream.
     *
     * @param chr1
     * @param chr2
     * @param reader
     * @param dis
     * @throws IOException
     */
    public MatrixZoomData(Chromosome chr1, Chromosome chr2, DatasetReader reader, LittleEndianInputStream dis) throws IOException {

        this.chr1 = chr1;
        this.chr2 = chr2;
        this.zoom = dis.readInt();

        if (reader.getVersion() >= 1) {
            dis.readInt();              // sum but we're not using this anymore
        }

        this.binSize = dis.readInt();
        this.blockBinCount = dis.readInt();
        this.blockColumnCount = dis.readInt();

        int nBlocks = dis.readInt();
        this.blockIndex = new HashMap<Integer, Preprocessor.IndexEntry>(nBlocks);

        for (int b = 0; b < nBlocks; b++) {
            int blockNumber = dis.readInt();
            long filePosition = dis.readLong();
            int blockSizeInBytes = dis.readInt();
            blockIndex.put(blockNumber, new Preprocessor.IndexEntry(filePosition, blockSizeInBytes));
        }

        blocks = new LinkedHashMap<Integer, Block>(nBlocks);
        this.reader = reader;

    }


    public int getBinSize() {
        return binSize;
    }

    public int getChr1() {
        return chr1.getIndex();
    }


    public int getChr2() {
        return chr2.getIndex();
    }

    public int getZoom() {
        return zoom;
    }

    public int getBlockColumnCount() {
        return blockColumnCount;
    }


    /**
     * Return the blocks overlapping the rectangular region specified.  The units are "bins"
     *
     * @param x1 leftmost position in "bins"
     * @param y1 top position in "bins"
     * @param x2 rightmost position in "bins"
     * @param y2 bottom position in "bins"
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

    public double[] getEigenvector() {
        return eigenvector;
    }

    public double[] computeEigenvector(DensityFunction df, int which) {

        if (pearsons == null) {
            pearsons = computePearsons(df);
        }
        int size = pearsons.getColumnDimension();
        eigenvector = new double[size];
        int numgood = 0;

        for (int i = 0; i < size; i++) {
            eigenvector[i] = Double.NaN;
            if (!isZeros(oe.getRow(i))) {
                eigenvector[i] = 1;
                numgood++;
            }
        }
        int[] cols = new int[numgood];
        numgood = 0;
        for (int i = 0; i < size; i++)
            if (!Double.isNaN(eigenvector[i]))
                cols[numgood++] = i;

        RealMatrix subMatrix = pearsons.getSubMatrix(cols, cols);

        if (which >= subMatrix.getColumnDimension() || which < 0)
            throw new NumberFormatException("Maximum eigenvector is " + subMatrix.getColumnDimension());


        RealVector rv = (new EigenDecompositionImpl(subMatrix, 0)).getEigenvector(which);
        double[] ev = rv.toArray();
        numgood = 0;
        for (int i = 0; i < size; i++) {
            if (Double.isNaN(eigenvector[i]))
                eigenvector[i] = 0;
            else
                eigenvector[i] = ev[numgood++];
        }

        return eigenvector;
    }

    public RealMatrix getPearsons() {
        return pearsons;
    }

    public RealMatrix computePearsons(DensityFunction df) {

        if (oe == null)
            oe = computeOE(df);


        pearsons = (new PearsonsCorrelation()).computeCorrelationMatrix(oe);

        PearsonsMinMax minMax = new PearsonsMinMax();
        pearsons.walkInOptimizedOrder(minMax);
        pearsonsMax = minMax.getMaxValue();
        pearsonsMin = minMax.getMinValue();
        return pearsons;
    }

    public double getPearsonsMin() {
        return pearsonsMin;
    }

    public double getPearsonsMax() {
        return pearsonsMax;
    }

    private RealMatrix readRealMatrix(String filename) throws IOException {
        LittleEndianInputStream is = null;
        RealMatrix rm = null;
        try {
            is = new LittleEndianInputStream(new BufferedInputStream(new FileInputStream(filename + this.zoom)));

            int rows = is.readInt();
            int cols = is.readInt();
            double[][] matrix = new double[rows][cols];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    matrix[i][j] = is.readDouble();
                }
            }
            rm = new Array2DRowRealMatrix(rows, cols);
            rm.setSubMatrix(matrix, 0, 0);
        } catch (IOException error) {
            System.err.println("IO error when saving Pearson's: " + error);
        } finally {
            if (is != null)
                is.close();
        }
        return rm;
    }

    private void outputRealMatrix(RealMatrix rm) throws IOException {
        LittleEndianOutputStream os = null;
        try {
            os = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream("pearsons" + this.zoom + ".bin")));

            int rows = rm.getRowDimension();
            int cols = rm.getColumnDimension();
            os.writeInt(rows);
            os.writeInt(cols);
            double[][] matrix = rm.getData();
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    os.writeDouble(matrix[i][j]);
                }
            }
        } catch (IOException error) {
            System.err.println("IO error when saving Pearson's: " + error);
        } finally {
            if (os != null)
                os.close();
        }
    }

    private boolean isZeros(double[] array) {
        for (double anArray : array)
            if (anArray != 0)
                return false;
        return true;
    }

    private double getVectorMean(RealVector vector) {
        double sum = 0;
        int count = 0;
        int size = vector.getDimension();
        for (int i = 0; i < size; i++) {
            if (!Double.isNaN(vector.getEntry(i))) {
                sum += vector.getEntry(i);
                count++;
            }
        }
        return sum / count;
    }

    public SparseRealMatrix computeOE(DensityFunction df) {

        if (chr1 != chr2) {
            throw new RuntimeException("Cannot yet compute Pearson's for different chromosomes");
        }

        int nBins = chr1.getSize() / binSize + 1;

        SparseRealMatrix rm = new OpenMapRealMatrix(nBins, nBins);

        List<Integer> blockNumbers = new ArrayList<Integer>(blockIndex.keySet());
        for (int blockNumber : blockNumbers) {
            Block b = readBlock(blockNumber);
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    int x = rec.getX();// * binSize;
                    int y = rec.getY();// * binSize;
                    int dist = Math.abs(x - y);
                    double expected = df.getDensity(chr1.getIndex(), dist);
                    double normCounts = (rec.getCounts() / expected);

                    rm.addToEntry(x, y, normCounts);
                    if (x != y) {
                        rm.addToEntry(y, x, normCounts);
                    }
                }
            }
        }
        int size = rm.getRowDimension();
        BitSet bitSet = new BitSet(size);
        double[] nans = new double[size];
        for (int i = 0; i < size; i++)
            nans[i] = Double.NaN;

        for (int i = 0; i < size; i++) {
            if (isZeros(rm.getRow(i))) {
                bitSet.set(i);
            }
        }

        for (int i = 0; i < size; i++) {
            if (bitSet.get(i)) {
                rm.setRow(i, nans);
                rm.setColumn(i, nans);
            }
        }

        for (int i = 0; i < size; i++) {
            RealVector v = rm.getRowVector(i);
            double m = getVectorMean(v);
            RealVector newV = v.mapSubtract(m);
            rm.setRowVector(i, newV);
        }
        PearsonsResetNan resetNan = new PearsonsResetNan();
        rm.walkInOptimizedOrder(resetNan);

        return rm;
    }

    /**
     * Compute scale parameters by from the first block of data
     *
     * @return scale parameters
     */
    public ScaleParameters computeScaleParameters() {
        double binSizeMB = binSize / 1000000.0;
        double binSizeMB2 = binSizeMB * binSizeMB;
        Block b = readBlock(0);
        if (b != null) {
            ContactRecord[] records = b.getContactRecords();
            double[] scores = new double[records.length];
            double sum = 0;
            for (int i = 0; i < scores.length; i++) {
                scores[i] = records[i].getCounts();
                sum += records[i].getCounts();
            }
            double percentile90 = StatUtils.percentile(scores, 90) / binSizeMB2;
            double mean = (sum / scores.length) / binSizeMB2;
            return new ScaleParameters(mean, percentile90);
        } else {
            return null;
        }

    }


    public void printDescription() {
        System.out.println("Chromosomes: " + chr1.getName() + " - " + chr2.getName());
        System.out.println("zoom: " + zoom);
        System.out.println("binSize (bp): " + binSize);
        System.out.println("blockBinCount (bins): " + blockBinCount);
        System.out.println("blockColumnCount (columns): " + blockColumnCount);

        System.out.println("Block size (bp): " + blockBinCount * binSize);
        System.out.println("");

    }
    // Dump the contents to standard out

    public void dump() {

        // Get the block index keys, and sort
        List<Integer> blockNumbers = new ArrayList<Integer>(blockIndex.keySet());
        Collections.sort(blockNumbers);

        System.out.println("# " + chr1.getName() + " - " + chr2.getName());

        for (int blockNumber : blockNumbers) {
            Block b = readBlock(blockNumber);
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    System.out.println(rec.getX() * binSize + "\t" + rec.getY() * binSize + "\t" + rec.getCounts());
                }
            }
        }
    }

    public void dumpOE(DensityFunction df) {
        oe = computeOE(df);
        int rows = oe.getRowDimension();
        int cols = oe.getColumnDimension();
        System.out.println(rows + " " + cols);
        double[][] matrix = oe.getData();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
        pearsons = (new PearsonsCorrelation()).computeCorrelationMatrix(oe);
        rows = pearsons.getRowDimension();
        cols = pearsons.getColumnDimension();
        System.out.println(rows + " " + cols);
        matrix = pearsons.getData();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();

    }

    private class PearsonsResetNan extends DefaultRealMatrixChangingVisitor {
        public double visit(int row, int column, double value) {
            if (Double.isNaN(value))
                return 0;
            return value;
        }
    }

    private class PearsonsMinMax extends DefaultRealMatrixPreservingVisitor {
        private double minValue = Double.MAX_VALUE;
        private double maxValue = Double.MIN_VALUE;

        public void visit(int row, int column, double value) {
            if (row != column) {
                if (value < minValue)
                    minValue = value;
                if (value > maxValue)
                    maxValue = value;
            }
        }

        public double getMinValue() {
            return minValue;
        }

        public double getMaxValue() {
            return maxValue;
        }

    }

}
