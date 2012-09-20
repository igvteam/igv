package org.broad.igv.hic.data;

import org.apache.commons.math.linear.*;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.matrix.BasicMatrix;
import org.broad.igv.hic.matrix.RealMatrixWrapper;
import org.broad.igv.hic.tools.Preprocessor;
import org.broad.igv.hic.track.HiCFixedGridAxis;
import org.broad.igv.hic.track.HiCGridAxis;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.collections.DoubleArrayList;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Aug 10, 2010
 */
public class MatrixZoomData {

    HiCGridAxis xGridAxis;
    HiCGridAxis yGridAxis;

    private Chromosome chr1;  // Redundant, but convenient
    private Chromosome chr2;  // Redundant, but convenient

    private int zoom;        // doesn't seem to be used
    private int fileBinSize;     // needed for fragment maps
    private int binSize;         // bin size in bp
    private int blockBinCount;   // block size in bins
    private int blockColumnCount;     // number of block columns

    // TODO -- isnt this a memory leak?  Should these be stored?
    private LinkedHashMap<Integer, Block> blocks;


    private Map<Integer, Preprocessor.IndexEntry> blockIndex;
    private DatasetReader reader;

    private BasicMatrix pearsons;
    private double[] eigenvector;
    private int[] nonCentromereColumns;

    public void setPearsons(BasicMatrix bm) {
        this.pearsons = bm;
    }

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

        fileBinSize = dis.readInt();

        this.blockBinCount = dis.readInt();
        this.blockColumnCount = dis.readInt();
        if (fileBinSize > 1) {
            binSize = fileBinSize;
        } else {
            binSize = (int) ((double) chr1.getLength() / (blockBinCount * blockColumnCount));
        }

        xGridAxis = new HiCFixedGridAxis(blockBinCount * blockColumnCount, binSize);
        yGridAxis = new HiCFixedGridAxis(blockBinCount * blockColumnCount, binSize);


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

        // If there's a pearson file available initialize it now
        String rootPath = FileUtils.getParent(reader.getPath());
        String folder = rootPath + "/" + chr1.getName();
        String file = "pearsons" + "_" + chr1.getName() + "_" + chr2.getName() + "_" + fileBinSize + ".bin";
        String fullPath = folder + "/" + file;
        if (FileUtils.resourceExists(fullPath)) {
            pearsons = ScratchPad.readPearsons(fullPath);
        }

        // If there's an eigenvector file load it
        String eigenFile = "eigen" + "_" + chr1.getName() + "_" + chr2.getName() + "_" + fileBinSize + ".wig";
        String fullEigenPath = folder + "/" + eigenFile;
        if (FileUtils.resourceExists(fullEigenPath)) {
            readEigenvector(fullEigenPath);
        }


    }


    public int getXBinCount() {
        return xGridAxis.getBinCount();
    }

    public int getYBinCount() {
        return yGridAxis.getBinCount();
    }

    public int getZoomMultiplier() {
        return binSize / 5000;
    }


    public Point getBinPosition(double genomePositionX, double genomePositionY) {

        int bX = xGridAxis.getBinNumberForGenomicPosition((int) genomePositionX);
        int bY = yGridAxis.getBinNumberForGenomicPosition((int) genomePositionY);
        return new Point(bX,  bY);

    }

    public Point getGenomePosition(int binX, int binY) {
        int gX = (xGridAxis.getGenomicStart(binX) + xGridAxis.getGenomicEnd(binX)) / 2;
        int gY = (yGridAxis.getGenomicStart(binY) + yGridAxis.getGenomicEnd(binY)) / 2;
        return new Point(gX, gY);

    }


    public int getMaxBinY() {
        return yGridAxis.getBinCount();
    }


    public int getMaxBinX() {
        return xGridAxis.getBinCount();
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

    private void readEigenvector(String fullPath) {

        if (FileUtils.resourceExists(fullPath)) {
            //Lots of assumptions made here about structure of wig file
            BufferedReader br = null;

            try {
                br = ParsingUtils.openBufferedReader(fullPath);
                String nextLine = br.readLine();  // The track line, ignored
                DoubleArrayList arrayList = new DoubleArrayList(10000);  // TODO -- can size this exactly
                while ((nextLine = br.readLine()) != null) {
                    if (nextLine.startsWith("track") || nextLine.startsWith("fixedStep") || nextLine.startsWith("#")) {
                        continue;
                    }
                    try {
                        arrayList.add(Double.parseDouble(nextLine));
                    } catch (NumberFormatException e) {
                        arrayList.add(Double.NaN);
                    }
                }
                eigenvector = arrayList.toArray();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } finally {
                if (br != null) try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

        }

    }


    public double[] computeEigenvector(DensityFunction df, int which) {
        //SparseRealMatrix oe = computeOE(df);
        if (pearsons == null) {
            pearsons = computePearsons(df);
        }

        RealMatrix subMatrix = null;
        if (pearsons instanceof RealMatrixWrapper) {

            subMatrix = ((RealMatrixWrapper) pearsons).getMatrix().getSubMatrix(nonCentromereColumns, nonCentromereColumns);

            if (which >= subMatrix.getColumnDimension() || which < 0)
                throw new NumberFormatException("Maximum eigenvector is " + subMatrix.getColumnDimension());


        } else {
            JOptionPane.showMessageDialog(null, "Eigenvector file not found");
            eigenvector = new double[0];
            return eigenvector;
            // TODO -- make submatrix from pearsons
            // Can't do this on precomputed matrix because we don't know the "nonCentromereColumns".
//            RealMatrix rm = new org.apache.commons.math.linear.BlockRealMatrix(pearsons.getRowDimension(), pearsons.getColumnDimension());
//            try {
//                for(int i=0; i<pearsons.getRowDimension(); i++) {
//                    for(int j=0; j<pearsons.getColumnDimension(); j++) {
//                        rm.setEntry(i, j, pearsons.getEntry(i, j));
//                    }
//                }
//            } catch (Exception e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
//
//            subMatrix = rm.getSubMatrix(nonCentromereColumns, nonCentromereColumns);

        }

        RealVector rv;
        rv = (new EigenDecompositionImpl(subMatrix, 0)).getEigenvector(which);

        double[] ev = rv.toArray();

        int size = pearsons.getColumnDimension();
        eigenvector = new double[size];
        int num = 0;
        for (int i = 0; i < size; i++) {
            if (i == nonCentromereColumns[num]) {
                eigenvector[i] = ev[num];
                num++;
            } else
                eigenvector[i] = 0;
        }
        return eigenvector;

    }

    public BasicMatrix getPearsons() {
        return pearsons;
    }

    public BasicMatrix computePearsons(DensityFunction df) {
        RealMatrix oe = computeOE(df);

        // below subtracts the empirical mean - necessary for mean-centered eigenvector
        int size = oe.getRowDimension();
        int num = 0;
        for (int i = 0; i < size; i++) {
            if (num < nonCentromereColumns.length && i == nonCentromereColumns[num]) {
                RealVector v = oe.getRowVector(i);
                double m = getVectorMean(v);
                RealVector newV = v.mapSubtract(m);
                oe.setRowVector(i, newV);
                num++;
            }
        }

        RealMatrix rm = (new PearsonsCorrelation()).computeCorrelationMatrix(oe);
        RealVector v = new ArrayRealVector(size);
        num = 0;
        for (int i = 0; i < size; i++) {
            if (num < nonCentromereColumns.length && i != nonCentromereColumns[num]) {
                rm.setRowVector(i, v);
                rm.setColumnVector(i, v);
            } else num++;
        }
        pearsons = new RealMatrixWrapper(rm);
        return pearsons;
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

        int nBins = blockBinCount * blockColumnCount;
        //int nBins = chr1.getLength() / binSize + 1;

        SparseRealMatrix rm = new OpenMapRealMatrix(nBins, nBins);

        List<Integer> blockNumbers = new ArrayList<Integer>(blockIndex.keySet());

        for (int blockNumber : blockNumbers) {
            Block b = readBlock(blockNumber);
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    int x = rec.getBinX();// * binSize;
                    int y = rec.getBinY();// * binSize;
                    int dist = Math.abs(x - y);
                    double expected = df.getDensity(chr1.getIndex(), dist);
                    double observed = rec.getCounts(); //df.getNormalizedCount(rec.getCounts(), chr1.getIndex(), x * binSize, chr2.getIndex(), y * binSize);
                    double normCounts = observed / expected;
                    rm.addToEntry(x, y, normCounts);
                    if (x != y) {
                        rm.addToEntry(y, x, normCounts);
                    }
                }
            }
        }
        int size = rm.getRowDimension();
        BitSet bitSet = new BitSet(size);

        for (int i = 0; i < size; i++) {
            if (isZeros(rm.getRow(i))) {
                bitSet.set(i);
            }
        }
        nonCentromereColumns = new int[size - bitSet.cardinality()];

        int num = 0;
        for (int i = 0; i < size; i++) {
            if (!bitSet.get(i)) {
                nonCentromereColumns[num++] = i;
            }
        }

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

    public void dump(LittleEndianOutputStream les) throws IOException {

        // Get the block index keys, and sort
        List<Integer> blockNumbers = new ArrayList<Integer>(blockIndex.keySet());
        Collections.sort(blockNumbers);

        if (les == null)
            System.out.println("# " + chr1.getName() + " - " + chr2.getName());


        for (int blockNumber : blockNumbers) {
            Block b = readBlock(blockNumber);
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    if (les == null)
                        System.out.println(rec.getBinX() * fileBinSize + "\t" + rec.getBinY() * fileBinSize + "\t" + rec.getCounts());
                    else {
                        les.writeInt(rec.getBinX());
                        les.writeInt(rec.getBinY());
                        les.writeInt(rec.getCounts());
                    }
                }
            }
        }
    }

    /**
     * Dump the O/E or Pearsons matrix to standard out in ascii format.
     *
     * @param df
     * @param isOE
     * @param les
     */
    public void dumpOE(DensityFunction df, boolean isOE, LittleEndianOutputStream les) throws IOException {

        if (isOE) {
            SparseRealMatrix oe = computeOE(df);
            int rows = oe.getRowDimension();
            int cols = oe.getColumnDimension();
            assert (rows == cols);
            if (les != null)
                les.writeInt(rows);
            else
                System.out.println(rows + " " + cols);
            int num = 0;
            for (int i = 0; i < rows; i++) {
                if (num >= nonCentromereColumns.length || i == nonCentromereColumns[num]) {
                    double[] row = oe.getRow(i);
                    int num2 = 0;
                    for (int j = 0; j < cols; j++) {
                        float output = Float.NaN;
                        if (num2 >= nonCentromereColumns.length || j == nonCentromereColumns[num2]) {
                            output = (float) row[j];
                            num2++;
                        }
                        if (les != null)
                            les.writeFloat(output);
                        else
                            System.out.print(output + " ");
                    }
                    num++;
                } else {
                    for (int j = 0; j < cols; j++) {
                        if (les != null)
                            les.writeFloat(Float.NaN);
                        else
                            System.out.print(Float.NaN + " ");
                    }
                }
                if (les == null)
                    System.out.println();
            }
            if (les == null)
                System.out.println();
        } else {

            RealMatrix rm = ((RealMatrixWrapper) computePearsons(df)).getMatrix();
            int rows = rm.getRowDimension();
            int cols = rm.getColumnDimension();
            if (les != null)
                les.writeInt(rows);
            else
                System.out.println(rows + " " + cols);
            double[][] matrix = rm.getData();
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    if (les != null)
                        les.writeFloat((float) matrix[i][j]);
                    else
                        System.out.print(matrix[i][j] + " ");
                }
                if (les == null)
                    System.out.println();
            }
            if (les == null)
                System.out.println();
        }
    }


}
