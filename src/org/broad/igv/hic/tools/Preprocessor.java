package org.broad.igv.hic.tools;

//import org.broad.igv.hic.MainWindow;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.data.*;
import org.broad.igv.util.CompressionUtils;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 * @since Aug 16, 2010
 */
public class Preprocessor {
    int nThreads = 0;

    private List<Chromosome> chromosomes;

    // Map of name -> index
    private Map<String, Integer> chromosomeOrdinals;

    private File outputFile;
    private LittleEndianOutputStream fos;

    private long masterIndexPosition;
    private Map<String, IndexEntry> matrixPositions;
    private Map<String, Long> blockIndexPositions;
    private Map<String, IndexEntry[]> blockIndexMap;

    private int countThreshold;
    private boolean diagonalsOnly;
    private boolean isNewVersion;
    private String fragmentFileName = null;
    private FragmentCalculation fragmentCalculation = null;
    private Set<String> includedChromosomes;

    public Preprocessor(File outputFile, List<Chromosome> chromosomes) {
        this.outputFile = outputFile;
        this.chromosomes = chromosomes;
        matrixPositions = new LinkedHashMap<String, IndexEntry>();
        blockIndexPositions = new LinkedHashMap<String, Long>();
        blockIndexMap = new LinkedHashMap<String, IndexEntry[]>();

        countThreshold = 0;
        diagonalsOnly = false;
        isNewVersion = false;
        chromosomeOrdinals = new Hashtable<String, Integer>();
        for (int i = 0; i < chromosomes.size(); i++) {
            chromosomeOrdinals.put(chromosomes.get(i).getName(), i);
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
    public void setNewVersion(boolean isNewVersion) {
        this.isNewVersion = isNewVersion;
    }

    public void setIncludedChromosomes(Set<String> includedChromosomes) {
        this.includedChromosomes = includedChromosomes;
    }

    public void setFragmentOption(String fragmentFileName) {
        this.fragmentFileName = fragmentFileName;
    }

    public void preprocess(final List<String> inputFileList) throws IOException {

        try {
            if (fragmentFileName != null) {
                fragmentCalculation = new FragmentCalculation(fragmentFileName, chromosomes);
            }
            else {
                System.out.println("WARNING: Not including fragment map");
            }

            System.out.println("Start preprocess");
            fos = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

            // Placeholder for master index position, replaced with actual position after all contents are written
            fos.writeLong(0l);

            int nChrs = chromosomes.size();
            fos.writeInt(nChrs);
            for (Chromosome chromosome : chromosomes) {
                fos.writeString(chromosome.getName());
                fos.writeInt(chromosome.getLength());
            }

            // Attribute dictionary -- only 1 attribute for now, version
            int nAttributes = 1;
            fos.writeInt(nAttributes);
            fos.writeString("Version");
            if (isNewVersion)          {
                fos.writeString("3");
            }
            else
                fos.writeString("2");

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
                    System.gc();
                }
            } // End of double loop through chromosomes


            masterIndexPosition = fos.getWrittenCount();
            writeMasterIndex(inputFileList);


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
     * Calculate observed/expected and write to a densities file that can be loaded later with
     * the Hi-C viewer.
     *
     * @param paths         Files to calculate densities on
     * @param buffer         Output stream for densities
     * @throws IOException
     */
    private void calculateDensities(List<String> paths, BufferedByteWriter buffer) throws IOException {
        DensityCalculation[] calcs = new DensityCalculation[Dataset.getNumberZooms()];
        for (int z = 0; z < Dataset.getNumberZooms(); z++) {
            if (fragmentCalculation == null && Dataset.getZoom(z) == 1)
                calcs[z] = null;
            else
                calcs[z] = new DensityCalculation(chromosomes, Dataset.getZoom(z), fragmentCalculation, isNewVersion);
        }
        if (isNewVersion) {
            System.out.println("Calculating coverage normalization");
            for (String path : paths) {
                PairIterator iter = (path.endsWith(".bin")) ?
                        new BinPairIterator(path) :
                        new AsciiPairIterator(path, chromosomeOrdinals);
                while (iter.hasNext()) {
                    AlignmentPair pair = iter.next();
                    for (int z = 0; z < Dataset.getNumberZooms(); z++) {
                        calcs[z].addToRow(pair);
                    }
                }
            }

            for (int z = 0; z < Dataset.getNumberZooms(); z++) {
                if (calcs[z] != null)
                    calcs[z].computeCoverageNormalization();
            }
        }
        System.out.println("Calculating distance normalization");
        for (String path : paths) {
            PairIterator iter = (path.endsWith(".bin")) ?
                    new BinPairIterator(path) :
                    new AsciiPairIterator(path, chromosomeOrdinals);
            while (iter.hasNext()) {
                AlignmentPair pair = iter.next();
                if (pair.getChr1() == pair.getChr2()) {
                    // Optionally filter on chromosome
                    if (includedChromosomes != null) {
                        String c1Name = chromosomes.get(pair.getChr1()).getName();
                        if (!(includedChromosomes.contains(c1Name))) {
                            continue;
                        }
                    }
                    int index = pair.getChr1();
                    for (int z = 0; z < Dataset.getNumberZooms(); z++) {
                        if (calcs[z] != null)
                            calcs[z].addDistance(index, pair.getPos1(), pair.getPos2());
                    }
                }
            }
        }
        for (int z = 0; z < Dataset.getNumberZooms(); z++) {
            if (calcs[z] != null)
                calcs[z].computeDensity();
        }
        System.out.println("Writing expected normalizations");
        outputDensities(calcs, buffer);

    }

    /**
     * Compute matrix for the given chromosome combination.  This results in full pass through the input files
     * for each chromosome combination.  This is done to save memory, at the expense of longer running times.
     *
     * @param inputFileList List of files to read
     * @param c1            Chromosome 1 -- always <= c2
     * @param c2            Chromosome 2
     * @return              Matrix with counts in each bin
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

        for (String file : inputFileList) {

            PairIterator iter = (file.endsWith(".bin")) ?
                    new BinPairIterator(file) :
                    new AsciiPairIterator(file, chromosomeOrdinals);

            while (iter.hasNext()) {

                AlignmentPair pair = iter.next();
                int pos1 = pair.getPos1();
                int pos2 = pair.getPos2();
                int chr1 = pair.getChr1();
                int chr2 = pair.getChr2();
                if (isWholeGenome) {
                    pos1 = getGenomicPosition(chr1, pos1);
                    pos2 = getGenomicPosition(chr2, pos2);
                    incrementCount(matrix, pos1, pos2);
                }
                else if ((c1 == chr1 && c2 == chr2) || (c1 == chr2 && c2 == chr1)) {
                    // we know c1 <= c2 and that's how the matrix is formed.
                    // pos1 goes with chr1 and pos2 goes with chr2
                    if (c1 == chr1)
                        incrementCount(matrix, pos1, pos2);
                    else // c1 == chr2
                        incrementCount(matrix, pos2, pos1);
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

    private static void incrementCount(MatrixPP matrix, int pos1,  int pos2) {
        matrix.incrementCount(pos1, pos2);
    }


    public void updateIndexPositions() throws IOException {
        RandomAccessFile raf = null;
        try {
            raf = new RandomAccessFile(outputFile, "rw");

            // Master index -- first entry in file (change later)
            raf.getChannel().position(0);
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

    public void writeMasterIndex(List<String> inputFileList) throws IOException {

        BufferedByteWriter buffer = new BufferedByteWriter();
        buffer.putInt(matrixPositions.size());
        for (Map.Entry<String, IndexEntry> entry : matrixPositions.entrySet()) {
            buffer.putString(entry.getKey());
            buffer.putLong(entry.getValue().position);
            buffer.putInt(entry.getValue().size);
        }
        calculateDensities(inputFileList, buffer);

        byte[] bytes = buffer.getBytes();
        fos.writeInt(bytes.length);
        fos.write(bytes);
    }

    private void outputDensities(DensityCalculation[] calcs, BufferedByteWriter buffer) throws IOException {
        int numCalcs = 0;
        for (DensityCalculation calc : calcs)
            if (calc != null)
                numCalcs++;

        buffer.putInt(numCalcs);
        for (DensityCalculation calc : calcs) {
            if (calc != null)
                calc.outputBinary(buffer, isNewVersion);
        }
    }


    public synchronized void writeMatrix(MatrixPP matrix) throws IOException {

        System.out.println("Start writing matrix: " + matrix.getKey());

        long position = fos.getWrittenCount();

        fos.writeInt(matrix.getChr1());
        fos.writeInt(matrix.getChr2());
        int numberZooms = 0;
        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            if (zd != null) {
                numberZooms++;
            }
        }
        fos.writeInt(numberZooms);
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
        return zd.getChr1() + "_" + zd.getChr2() + "_" + zd.getZoom();
    }

    private void writeZoomHeader(MatrixZoomDataPP zd) throws IOException {

        int numberOfBlocks = zd.getBlocks().size();

        fos.writeInt(zd.getZoom());
        fos.writeInt(zd.getSum());
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
                buffer.putInt(rec.getCounts());
            }
        }

        byte[] bytes = buffer.getBytes();
        byte[] compressedBytes = CompressionUtils.compress(bytes);
        fos.write(compressedBytes);

    }


    static public class BufferedByteWriter {

        ByteArrayOutputStream buffer;
        LittleEndianOutputStream dos;

        public BufferedByteWriter() {
            this(8192);
        }


        public BufferedByteWriter(int size) {
            if (size <= 0) {
                throw new IllegalArgumentException("Buffer size <= 0");
            }
            buffer = new ByteArrayOutputStream(size);
            dos = new LittleEndianOutputStream(buffer);
        }

        public byte[] getBytes() {
            return buffer.toByteArray();
        }

        private void put(byte[] b) throws IOException {
            dos.write(b);
        }

        private void put(byte b) throws IOException {
            dos.write(b);
        }

        private void putShort(short v) throws IOException {

            dos.writeShort(v);
        }

        public void putInt(int v) throws IOException {
            dos.writeInt(v);
        }

        public void putDouble(double v) throws IOException {
            dos.writeDouble(v);
        }

        public void putLong(long v) throws IOException {
            dos.writeLong(v);
        }

        public void putString(String string) throws IOException {
            dos.writeString(string);
        }

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

        private int chr1;
        private int chr2;
        private MatrixZoomDataPP[] zoomData;


        /**
         * Constructor for creating a matrix and initializing zoomed data at predefined resolution scales.  This
         * constructor is used when parsing alignment files.
         *
         * @param chr1   Chromosome 1
         * @param chr2   Chromosome 2
         */
        MatrixPP(int chr1, int chr2) {
            this.chr1 = chr1;
            this.chr2 = chr2;
            zoomData = new MatrixZoomDataPP[Dataset.getNumberZooms()];
            for (int zoom = 0; zoom < Dataset.getNumberZooms(); zoom++) {
                int binSize = Dataset.getZoom(zoom);
                Chromosome chrom1 = chromosomes.get(chr1);
                Chromosome chrom2 = chromosomes.get(chr2);
                int nBins;

                if (binSize == 1 && fragmentCalculation != null) {
                    // Size block (submatrices) to be ~500 bins wide.
                    nBins = Math.max(fragmentCalculation.getNumberFragments(chrom1), fragmentCalculation.getNumberFragments(chrom2));
                    int nColumns = Math.max(1, nBins / 500);
                    zoomData[zoom] = new MatrixZoomDataPP(chrom1, chrom2, binSize, nColumns, zoom);

                }
                else  if (binSize > 1) {       // this might break in the meantime, before frags are completely implemented.
                    // Size block (submatrices) to be ~500 bins wide.
                    int len = Math.max(chrom1.getLength(), chrom2.getLength());
                    nBins = len / binSize;   // Size of chrom in bins
                    int nColumns = Math.max(1, nBins / 500);
                    zoomData[zoom] = new MatrixZoomDataPP(chrom1, chrom2, binSize, nColumns, zoom);

                }

            }
        }

        /**
         * Constructor for creating a matrix with a single zoom level at a specified bin size.  This is provided
         * primarily for constructing a whole-genome view.
         *
         * @param chr1   Chromosome 1
         * @param chr2   Chromosome 2
         * @param binSize   Bin size
         */
        MatrixPP(int chr1, int chr2, int binSize) {
            this.chr1 = chr1;
            this.chr2 = chr2;
            zoomData = new MatrixZoomDataPP[1];
            int nBlocks = 1;
            zoomData[0] = new MatrixZoomDataPP(chromosomes.get(chr1), chromosomes.get(chr2), binSize, nBlocks, 0);

        }


        String generateKey(int chr1, int chr2) {
            return "" + chr1 + "_" + chr2;
        }

        String getKey() {
            return generateKey(chr1, chr2);
        }


        void incrementCount(int pos1, int pos2) {

            for (MatrixZoomDataPP aZoomData : zoomData) {
                if (aZoomData != null)   // fragment level could be null
                    aZoomData.incrementCount(pos1, pos2);
            }

        }

        void parsingComplete() {
            for (MatrixZoomDataPP zd : zoomData) {
                if (zd != null) // fragment level could be null
                    zd.parsingComplete();
            }
        }

        int getChr1() {
            return chr1;
        }

        int getChr2() {
            return chr2;
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

        private int sum;
        private int zoom;
        private int binSize;         // bin size in bp
        private int blockBinCount;   // block size in bins
        private int blockColumnCount;     // number of block columns

        private LinkedHashMap<Integer, Block> blocks;
        private Map<Integer, Preprocessor.IndexEntry> blockIndex;

        int getSum() {
            return sum;
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
         * Representation of MatrixZoomData used for preprocessing
         *
         * @param chr1             index of first chromosome  (x-axis)
         * @param chr2             index of second chromosome
         * @param binSize          size of each grid bin in bp
         * @param blockColumnCount number of block columns
         * @param zoom             integer zoom (resolution) level index.  TODO Is this needed?
         */
        MatrixZoomDataPP(Chromosome chr1, Chromosome chr2, int binSize, int blockColumnCount, int zoom) {

            this.sum = 0;
            this.chr1 = chr1;
            this.chr2 = chr2;
            this.binSize = binSize;
            this.blockColumnCount = blockColumnCount;
            this.zoom = zoom;

            int nBinsX;
            if (binSize == 1) {
                nBinsX = fragmentCalculation.getNumberFragments(chr1);
            }
            else
                nBinsX = chr1.getLength() / binSize + 1;
            blockBinCount = nBinsX / blockColumnCount + 1;
            blocks = new LinkedHashMap<Integer, Block>(blockColumnCount * blockColumnCount);
        }


        public void incrementCount(int pos1, int pos2) {

            sum++;

            int xBin;
            int yBin;

            if (binSize == 1) {

                xBin = fragmentCalculation.getBin(chr1, pos1);
                yBin = fragmentCalculation.getBin(chr2, pos2);
//                if (chr1.getIndex() == 21 && xBin > 9000)
//                    System.out.println(chr1 + " " + pos1 + " " + xBin + " " + chr2 + " " + pos2 + " " + yBin);
//                if (chr2.getIndex() == 21 && yBin > 9000)
//                    System.out.println(chr1 + " " + pos1 + " " + xBin + " " + chr2 + " " + pos2 + " " + yBin);
            }
            else {
                xBin = pos1 / getBinSize();
                yBin = pos2 / getBinSize();
            }
            if (chr1.equals(chr2)) {
                int b1 = Math.min(xBin, yBin);
                int b2 = Math.max(xBin, yBin);
                xBin = b1;
                yBin = b2;

                if (b1 != b2) {
                    sum++;  // <= count for mirror cell.
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
            block.incrementCount(xBin, yBin);
        }

        void parsingComplete() {
            for (Block b : blocks.values()) {
                b.parsingComplete();
            }
        }


    }
}
