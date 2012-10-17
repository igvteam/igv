package org.broad.igv.hic.tools;

//import org.broad.igv.hic.MainWindow;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.data.*;
import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.util.CompressionUtils;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 * @since Aug 16, 2010
 */
public class Preprocessor {

    // Base-pair resultions
    public static final int[] bpBinSizes = {2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000};
    public static final String[] bpResLabels = {"2.5 MB", "1 MB", "500 KB", "250 KB", "100 KB", "50 KB", "25 KB", "10 KB", "5 KB"};

    // Fragment resolutions
    public static final int[] fragBinSizes = {500, 200, 100, 50, 20, 5, 2, 1};
    public static final String[] fragResLabels = {"500f", "200f", "100f", "50f", "20f", "5f", "2f", "1f"};

    int nThreads = 0;

    private List<Chromosome> chromosomes;

    // Map of name -> index
    private Map<String, Integer> chromosomeIndexes;

    private File outputFile;
    private LittleEndianOutputStream fos;

    private long masterIndexPosition;
    private Map<String, IndexEntry> matrixPositions;
    private Map<String, Long> blockIndexPositions;
    private Map<String, IndexEntry[]> blockIndexMap;

    private int countThreshold = 0;
    private boolean diagonalsOnly = false;
    private String fragmentFileName = null;
    private FragmentCalculation fragmentCalculation = null;
    private Set<String> includedChromosomes;
    private String genomeId;

    /**
     * The position of the field containing the masterIndex position
     */
    private long masterIndexPositionPosition;

    private Map<String, ExpectedValueCalculation> expectedValueCalculations;

    public Preprocessor(File outputFile, String genomeId, List<Chromosome> chromosomes) {
        this.genomeId = genomeId;
        this.outputFile = outputFile;
        this.matrixPositions = new LinkedHashMap<String, IndexEntry>();
        this.blockIndexPositions = new LinkedHashMap<String, Long>();
        this.blockIndexMap = new LinkedHashMap<String, IndexEntry[]>();

        this.chromosomes = chromosomes;
        chromosomeIndexes = new Hashtable<String, Integer>();
        for (int i = 0; i < chromosomes.size(); i++) {
            chromosomeIndexes.put(chromosomes.get(i).getName(), i);
        }

    }

    public void setNumberOfThreads(int n) {
        this.nThreads = n;
    }

    public void setCountThreshold(int countThreshold) {
        this.countThreshold = countThreshold;
    }

    public void setDiagonalsOnly(boolean diagonalsOnly) {
        this.diagonalsOnly = diagonalsOnly;
    }

    public void setIncludedChromosomes(Set<String> includedChromosomes) {
        this.includedChromosomes = includedChromosomes;
    }

    public void setFragmentFile(String fragmentFileName) {
        this.fragmentFileName = fragmentFileName;
    }

    public void preprocess(final List<String> inputFileList) throws IOException {

        try {
            if (fragmentFileName != null) {
                fragmentCalculation = new FragmentCalculation(fragmentFileName, chromosomes);
            } else {
                System.out.println("WARNING: Not including fragment map");
            }

            expectedValueCalculations = new LinkedHashMap<String, ExpectedValueCalculation>();
            for (int z = 0; z < bpBinSizes.length; z++) {
                ExpectedValueCalculation calc = new ExpectedValueCalculation(chromosomes, bpBinSizes[z], null);
                String key = "BP_" + bpBinSizes[z];
                expectedValueCalculations.put(key, calc);
            }
            if (fragmentCalculation != null) {
                for (int z = 0; z < fragBinSizes.length; z++) {
                    ExpectedValueCalculation calc = new ExpectedValueCalculation(chromosomes, fragBinSizes[z], fragmentCalculation);
                    String key = "FRAG_" + fragBinSizes[z];
                    expectedValueCalculations.put(key, calc);
                }
            }


            System.out.println("Start preprocess");
            fos = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

            // Magic number
            byte[] magicBytes = "HIC".getBytes();
            fos.write(magicBytes[0]);
            fos.write(magicBytes[1]);
            fos.write(magicBytes[2]);
            fos.write(0);

            // Version
            fos.writeInt(4);

            // Placeholder for master index position, replaced with actual position after all contents are written
            masterIndexPositionPosition = fos.getWrittenCount();
            fos.writeLong(0l);


            // Genome ID
            fos.writeString(genomeId);

            // Sequence dictionary
            int nChrs = chromosomes.size();
            fos.writeInt(nChrs);
            for (Chromosome chromosome : chromosomes) {
                fos.writeString(chromosome.getName());
                fos.writeInt(chromosome.getLength());
            }

            //BP resolution levels
            int nBpRes = bpBinSizes.length;
            fos.writeInt(nBpRes);
            for (int i = 0; i < nBpRes; i++) {
                fos.writeInt(bpBinSizes[i]);
            }

            //fragment resolutions
            int nFragRes = fragmentCalculation == null ? 0 : fragBinSizes.length;
            fos.writeInt(nFragRes);
            for (int i = 0; i < nFragRes; i++) {
                fos.writeInt(fragBinSizes[i]);
            }

            // fragment sites
            if (nFragRes > 0) {
                for (Chromosome chromosome : chromosomes) {
                    int[] sites = fragmentCalculation.getSites(chromosome.getName());
                    int nSites = sites == null ? 0 : sites.length;
                    fos.writeInt(nSites);
                    for (int i = 0; i < nSites; i++) {
                        fos.writeInt(sites[i]);
                    }
                }
            }

            //hemi-frag levels (reserved for future use)
            int nHemiFrags = 0;
            fos.writeInt(nHemiFrags);


            // Attribute dictionary
            int nAttributes = 0;
            fos.writeInt(nAttributes);
            //fos.writeString("theKey");
            //fos.writeString("theValue");
            //... repeat for each attribute

            // Compute matrices.  Note that c2 is always >= c1
            for (int c1 = 0; c1 < nChrs; c1++) {
                for (int c2 = c1; c2 < nChrs; c2++) {

                    // Index zero is whole genome
                    if ((c1 == 0 && c2 != 0) || (c2 == 0 && c1 != 0)) continue;

                    if (diagonalsOnly && c1 != c2) continue;

                    // Optionally filter on chromosome
                    if (includedChromosomes != null && c1 != 0) {
                        String c1Name = chromosomes.get(c1).getName();
                        String c2Name = chromosomes.get(c2).getName();
                        if (!(includedChromosomes.contains(c1Name) || includedChromosomes.contains(c2Name))) {
                            continue;
                        }
                    }

                    MatrixPP matrix = computeMatrix(inputFileList, c1, c2);
                    if (matrix != null) {
                        writeMatrix(matrix);
                    }
                }
            } // End of double loop through chromosomes


            masterIndexPosition = fos.getWrittenCount();
            writeFooter();


        } finally {
            if (fos != null)
                fos.close();
        }

        updateIndexPositions();
    }

    private void processThreads(List<Thread> threads, List<MatrixPP> matrices) throws IOException {
        for (Thread t : threads) {
            t.start();
        }
        for (Thread t : threads) {
            try {
                t.join();
            } catch (InterruptedException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        for (MatrixPP matrix : matrices) {
            writeMatrix(matrix);
        }
    }


    /**
     * Compute matrix for the given chromosome combination.  This results in full pass through the input files
     * for each chromosome combination.  This is done to save memory, at the expense of longer running times.
     *
     * @param inputFileList List of files to read
     * @param c1            Chromosome 1 -- always <= c2
     * @param c2            Chromosome 2
     * @return Matrix with counts in each bin
     * @throws IOException
     */
    public MatrixPP computeMatrix(List<String> inputFileList, int c1, int c2) throws IOException {

        boolean isWholeGenome = (c1 == 0 && c2 == 0);

        MatrixPP matrix;
        // NOTE: always true that c1 <= c2
        if (isWholeGenome) {
            int genomeLength = chromosomes.get(0).getLength();  // <= whole genome in KB
            int binSize = genomeLength / 500;
            matrix = new MatrixPP(c1, c2, binSize);
        } else {
            matrix = new MatrixPP(c1, c2);
        }

        // TODO -- in the future this value will be read from the file, might not be 1
        float score = 1.0f;

        for (String file : inputFileList) {

            PairIterator iter = (file.endsWith(".bin")) ?
                    new BinPairIterator(file, chromosomeIndexes) :
                    new AsciiPairIterator(file, chromosomeIndexes);

            while (iter.hasNext()) {

                AlignmentPair pair = iter.next();
                int pos1 = pair.getPos1();
                int pos2 = pair.getPos2();
                int chr1 = pair.getChr1();
                int chr2 = pair.getChr2();
                if (isWholeGenome) {
                    pos1 = getGenomicPosition(chr1, pos1);
                    pos2 = getGenomicPosition(chr2, pos2);
                    matrix.incrementCount(pos1, pos2, score);
                } else if ((c1 == chr1 && c2 == chr2) || (c1 == chr2 && c2 == chr1)) {
                    // we know c1 <= c2 and that's how the matrix is formed.
                    // pos1 goes with chr1 and pos2 goes with chr2
                    if (c1 == chr1) {
                        matrix.incrementCount(pos1, pos2, score);
                    } else {// c1 == chr2
                        matrix.incrementCount(pos2, pos1, score);
                    }

                }

            }

            iter.close();

        }


        matrix.parsingComplete();
        return matrix;
    }

    private int getGenomicPosition(int chr, int pos) {
        long len = 0;
        for (int i = 1; i < chr; i++) {
            len += chromosomes.get(i).getLength();
        }
        len += pos;

        return (int) (len / 1000);

    }

    public void updateIndexPositions() throws IOException {
        RandomAccessFile raf = null;
        try {
            raf = new RandomAccessFile(outputFile, "rw");

            // Master index
            raf.getChannel().position(masterIndexPositionPosition);
            BufferedByteWriter buffer = new BufferedByteWriter();
            buffer.putLong(masterIndexPosition);
            raf.write(buffer.getBytes());

            // Block indices
            for (String key : blockIndexPositions.keySet()) {
                long pos = blockIndexPositions.get(key);
                IndexEntry[] blockIndex = blockIndexMap.get(key);

                if (blockIndex == null) {
                    System.err.println("Missing block index for: " + key);
                } else {
                    raf.getChannel().position(pos);

                    // Write as little endian
                    buffer = new BufferedByteWriter();
                    for (IndexEntry aBlockIndex : blockIndex) {
                        buffer.putInt(aBlockIndex.id);
                        buffer.putLong(aBlockIndex.position);
                        buffer.putInt(aBlockIndex.size);
                    }
                    raf.write(buffer.getBytes());
                }
            }
        } finally {
            if (raf != null) raf.close();
        }
    }

    public void writeFooter() throws IOException {

        BufferedByteWriter buffer = new BufferedByteWriter();
        buffer.putInt(matrixPositions.size());
        for (Map.Entry<String, IndexEntry> entry : matrixPositions.entrySet()) {
            buffer.putNullTerminatedString(entry.getKey());
            buffer.putLong(entry.getValue().position);
            buffer.putInt(entry.getValue().size);
        }


        buffer.putInt(expectedValueCalculations.size());
        for (Map.Entry<String, ExpectedValueCalculation> entry : expectedValueCalculations.entrySet()) {
            String key = entry.getKey();
            ExpectedValueCalculation ev = entry.getValue();

            ev.computeDensity();

            System.out.println(key + "  Norm factors");

            // Map of chromosome index -> normalization factor
            Map<Integer, Double> normalizationFactors = ev.getNormalizationFactors();
            buffer.putInt(normalizationFactors.size());
            for (Map.Entry<Integer, Double> normFactor : normalizationFactors.entrySet()) {
                buffer.putInt(normFactor.getKey());
                buffer.putDouble(normFactor.getValue());
                System.out.println(normFactor.getKey() + "  " + normFactor.getValue());
            }

            // The density values
            double[] expectedValues = ev.getDensityAvg();
            buffer.putNullTerminatedString(key);
            buffer.putInt(expectedValues.length);
            for (int i = 0; i < expectedValues.length; i++) {
                buffer.putDouble(expectedValues[i]);
            }
        }

        byte[] bytes = buffer.getBytes();
        fos.writeInt(bytes.length);
        fos.write(bytes);
    }

    public synchronized void writeMatrix(MatrixPP matrix) throws IOException {

        System.out.println("Start writing matrix: " + matrix.getKey());

        long position = fos.getWrittenCount();

        fos.writeInt(matrix.getChr1Idx());
        fos.writeInt(matrix.getChr2Idx());
        int numResolutions = 0;

        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            if (zd != null) {
                numResolutions++;
            }
        }
        fos.writeInt(numResolutions);

        //fos.writeInt(matrix.getZoomData().length);
        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            if (zd != null)
                writeZoomHeader(zd);
        }

        int size = (int) (fos.getWrittenCount() - position);
        matrixPositions.put(matrix.getKey(), new IndexEntry(position, size));

        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            if (zd != null) {
                IndexEntry[] blockIndex = writeZoomData(zd);
                final String blockKey = getBlockKey(zd);
                blockIndexMap.put(blockKey, blockIndex);
            }
        }

        System.out.println("Done writing matrix: " + matrix.getKey());
    }

    private String getBlockKey(MatrixZoomDataPP zd) {
        return zd.getChr1() + "_" + zd.getChr2() + "_" + zd.getUnit() + "_" + zd.getBinSize();
    }

    private void writeZoomHeader(MatrixZoomDataPP zd) throws IOException {

        int numberOfBlocks = zd.getBlocks().size();

        fos.writeString(zd.getUnit());  // Unit, ether "BP" or "FRAG"
        fos.writeInt(zd.getZoom());     // zoom index,  lowest res is zero
        fos.writeFloat((float) zd.getSum());      // sum
        fos.writeFloat((float) zd.getAverage());
        fos.writeFloat((float) zd.getStdDev());
        fos.writeFloat((float) zd.getPercent95());
        fos.writeInt(zd.getBinSize());
        fos.writeInt(zd.getBlockBinCount());
        fos.writeInt(zd.getBlockColumnCount());
        fos.writeInt(numberOfBlocks);
        blockIndexPositions.put(getBlockKey(zd), fos.getWrittenCount());

        // Placeholder for block index
        for (int i = 0; i < numberOfBlocks; i++) {
            fos.writeInt(0);
            fos.writeLong(0l);
            fos.writeInt(0);
        }

    }

    private IndexEntry[] writeZoomData(MatrixZoomDataPP zd) throws IOException {

        final Map<Integer, Block> blocks = zd.getBlocks();

        IndexEntry[] indexEntries = new IndexEntry[blocks.size()];
        int i = 0;
        for (Map.Entry<Integer, Block> entry : blocks.entrySet()) {

            int blockNumber = entry.getKey();
            Block block = entry.getValue();

            long position = fos.getWrittenCount();
            writeContactRecords(block);
            int size = (int) (fos.getWrittenCount() - position);

            indexEntries[i] = new IndexEntry(blockNumber, position, size);
            i++;
        }
        return indexEntries;

    }

    /**
     * Note -- compressed
     *
     * @param block Block to write
     * @throws IOException
     */
    private void writeContactRecords(Block block) throws IOException {

        final Collection<ContactRecord> records = block.getContractRecordValues();//   getContactRecords();


        // Count records first
        int nRecords;
        if (countThreshold > 0) {
            nRecords = 0;
            for (ContactRecord rec : records) {
                if (rec.getCounts() >= countThreshold) {
                    nRecords++;
                }
            }
        } else {
            nRecords = records.size();
        }

        BufferedByteWriter buffer = new BufferedByteWriter(nRecords * 12);
        buffer.putInt(nRecords);
        for (ContactRecord rec : records) {
            if (rec.getCounts() >= countThreshold) {
                buffer.putInt(rec.getBinX());
                buffer.putInt(rec.getBinY());
                buffer.putFloat(rec.getCounts());
            }
        }

        byte[] bytes = buffer.getBytes();
        byte[] compressedBytes = CompressionUtils.compress(bytes);
        fos.write(compressedBytes);

    }


    public static class IndexEntry {
        int id;
        public long position;
        public int size;

        IndexEntry(int id, long position, int size) {
            this.id = id;
            this.position = position;
            this.size = size;
        }

        public IndexEntry(long position, int size) {
            this.position = position;
            this.size = size;
        }
    }

    /**
     * @author jrobinso
     * @since Aug 12, 2010
     */
    class MatrixPP {

        private int chr1Idx;
        private int chr2Idx;
        private MatrixZoomDataPP[] zoomData;


        /**
         * Constructor for creating a matrix and initializing zoomed data at predefined resolution scales.  This
         * constructor is used when parsing alignment files.
         *
         * @param chr1Idx Chromosome 1
         * @param chr2Idx Chromosome 2
         */
        MatrixPP(int chr1Idx, int chr2Idx) {
            this.chr1Idx = chr1Idx;
            this.chr2Idx = chr2Idx;

            int nResolutions = bpBinSizes.length;
            if (fragmentCalculation != null) {
                nResolutions += fragBinSizes.length;
            }

            zoomData = new MatrixZoomDataPP[nResolutions];

            int zoom = 0; //
            for (int idx = 0; idx < bpBinSizes.length; idx++) {
                int binSize = bpBinSizes[zoom];
                Chromosome chrom1 = chromosomes.get(chr1Idx);
                Chromosome chrom2 = chromosomes.get(chr2Idx);

                // Size block (submatrices) to be ~500 bins wide.
                int len = Math.max(chrom1.getLength(), chrom2.getLength());
                int nBins = len / binSize;   // Size of chrom in bins
                int nColumns = Math.max(1, nBins / 500);
                zoomData[idx] = new MatrixZoomDataPP(chrom1, chrom2, binSize, nColumns, zoom, false);
                zoom++;

            }

            if (fragmentCalculation != null) {
                Chromosome chrom1 = chromosomes.get(chr1Idx);
                Chromosome chrom2 = chromosomes.get(chr2Idx);
                int nFragBins1 = Math.max(fragmentCalculation.getNumberFragments(chrom1.getName()),
                        fragmentCalculation.getNumberFragments(chrom2.getName()));

                zoom = 0;
                for (int idx = bpBinSizes.length; idx < nResolutions; idx++) {
                    int binSize = fragBinSizes[zoom];
                    int nBins = nFragBins1 / binSize;
                    int nColumns = Math.max(1, nBins / 500);
                    zoomData[idx] = new MatrixZoomDataPP(chrom1, chrom2, binSize, nColumns, zoom, true);
                    zoom++;
                }
            }
        }

        /**
         * Constructor for creating a matrix with a single zoom level at a specified bin size.  This is provided
         * primarily for constructing a whole-genome view.
         *
         * @param chr1Idx    Chromosome 1
         * @param chr2Idx    Chromosome 2
         * @param binSize Bin size
         */
        MatrixPP(int chr1Idx, int chr2Idx, int binSize) {
            this.chr1Idx = chr1Idx;
            this.chr2Idx = chr2Idx;
            zoomData = new MatrixZoomDataPP[1];
            int nBlocks = 1;
            zoomData[0] = new MatrixZoomDataPP(chromosomes.get(chr1Idx), chromosomes.get(chr2Idx), binSize, nBlocks, 0, false);

        }


        String getKey() {
            return "" + chr1Idx + "_" + chr2Idx;
        }


        void incrementCount(int pos1, int pos2, float score) {

            for (MatrixZoomDataPP aZoomData : zoomData) {
                aZoomData.incrementCount(pos1, pos2, score);
            }
        }

        void parsingComplete() {
            for (MatrixZoomDataPP zd : zoomData) {
                if (zd != null) // fragment level could be null
                    zd.parsingComplete();
            }
        }

        int getChr1Idx() {
            return chr1Idx;
        }

        int getChr2Idx() {
            return chr2Idx;
        }

        MatrixZoomDataPP[] getZoomData() {
            return zoomData;
        }

    }


    /**
     * @author jrobinso
     * @since Aug 10, 2010
     */
    class MatrixZoomDataPP {

        private Chromosome chr1;  // Redundant, but convenient    BinDatasetReader
        private Chromosome chr2;  // Redundant, but convenient

        private double sum;
        private double average;    // Number of occupied cells;
        private double stdDev;
        private double percent95;

        private int zoom;
        private int binSize;         // bin size in bp
        private int blockBinCount;   // block size in bins
        private int blockColumnCount;     // number of block columns

        boolean isFrag;

        private LinkedHashMap<Integer, Block> blocks;
        private Map<Integer, Preprocessor.IndexEntry> blockIndex;


        /**
         * Representation of MatrixZoomData used for preprocessing
         *
         * @param chr1             index of first chromosome  (x-axis)
         * @param chr2             index of second chromosome
         * @param binSize          size of each grid bin in bp
         * @param blockColumnCount number of block columns
         * @param zoom             integer zoom (resolution) level index.  TODO Is this needed?
         */
        MatrixZoomDataPP(Chromosome chr1, Chromosome chr2, int binSize, int blockColumnCount, int zoom, boolean isFrag) {

            this.sum = 0;
            this.chr1 = chr1;
            this.chr2 = chr2;
            this.binSize = binSize;
            this.blockColumnCount = blockColumnCount;
            this.zoom = zoom;
            this.isFrag = isFrag;

            // Get length in proper units
            int len = isFrag ? fragmentCalculation.getNumberFragments(chr1.getName()) : chr1.getLength();

            int nBinsX = len / binSize + 1;

            blockBinCount = nBinsX / blockColumnCount + 1;
            blocks = new LinkedHashMap<Integer, Block>(blockColumnCount * blockColumnCount);
        }

        String getUnit() {
            return isFrag ? "FRAG" : "BP";
        }

        double getSum() {
            return sum;
        }

        double getAverage() {
            return average;
        }

        public double getStdDev() {
            return stdDev;
        }

        public double getPercent95() {
            return percent95;
        }

        private void computeStats() {
            int cellCount = 0;
            average = 0;
            double sum = 0;
            for (Block b : blocks.values()) {
                for (ContactRecord cr : b.getContractRecordValues()) {
                    sum += cr.getCounts();
                    cellCount++;
                }
            }
            average = sum / cellCount;

            double sumSq = 0;
            int idx = 0;
            double[] values = new double[cellCount];
            for (Block b : blocks.values()) {
                for (ContactRecord cr : b.getContractRecordValues()) {
                    double diff = cr.getCounts() - average;
                    sumSq += (diff * diff);
                    values[idx++] = cr.getCounts();
                }
            }
            stdDev = Math.sqrt(sumSq / cellCount);
            percent95 = StatUtils.percentile(values, 90);
        }


        int getBinSize() {
            return binSize;
        }


        Chromosome getChr1() {
            return chr1;
        }


        Chromosome getChr2() {
            return chr2;
        }

        int getZoom() {
            return zoom;
        }

        int getBlockBinCount() {
            return blockBinCount;
        }

        int getBlockColumnCount() {
            return blockColumnCount;
        }

        Map<Integer, Block> getBlocks() {
            return blocks;
        }

        /**
         * Increment the count for the bin represented by the GENOMIC position (pos1, pos2)
         *
         * @param gPos1
         * @param gPos2
         */
        public void incrementCount(int gPos1, int gPos2, float score) {

            sum += score;
            // Convert to proper units,  fragments or base-pairs

            int pos1 = isFrag ? fragmentCalculation.getBin(chr1.getName(), gPos1) : gPos1;
            int pos2 = isFrag ? fragmentCalculation.getBin(chr2.getName(), gPos2) : gPos2;
            if (pos1 < 0 || pos2 < 0) return;

            int xBin = pos1 / binSize;
            int yBin = pos2 / binSize;

            if (chr1.equals(chr2)) {
                int b1 = Math.min(xBin, yBin);
                int b2 = Math.max(xBin, yBin);
                xBin = b1;
                yBin = b2;

                if (b1 != b2) {
                    sum++;  // <= count for mirror cell.
                }

                String evKey = (isFrag ? "FRAG_" : "BP_") + binSize;
                ExpectedValueCalculation ev = expectedValueCalculations.get(evKey);
                if (ev != null) {
                    ev.addDistance(chr1.getIndex(), xBin, yBin);
                }
            }

            // compute block number (fist block is zero)
            int blockCol = xBin / getBlockBinCount();
            int blockRow = yBin / getBlockBinCount();
            int blockNumber = getBlockColumnCount() * blockRow + blockCol;

            Block block = blocks.get(blockNumber);
            if (block == null) {
                block = new Block(blockNumber);
                blocks.put(blockNumber, block);
            }
            block.incrementCount(xBin, yBin, score);
        }

        void parsingComplete() {
            for (Block b : blocks.values()) {
                b.parsingComplete();

            }
            computeStats();
        }


    }
}
