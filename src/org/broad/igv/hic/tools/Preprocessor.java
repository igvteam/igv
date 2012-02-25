package org.broad.igv.hic.tools;

//import org.broad.igv.hic.MainWindow;
import org.broad.igv.hic.HiCGlobals;
import org.broad.igv.hic.data.*;
import org.broad.igv.util.CompressionUtils;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 * @date Aug 16, 2010
 */
public class Preprocessor {

    List<Chromosome> chromosomes;

    // Map of name -> index
    Map<String, Integer> chromosomeOrdinals;

    File outputFile;
    LittleEndianOutputStream fos;
    long bytesWritten = 0;
    long totalCount = 0;

    long masterIndexPosition;
    Map<String, IndexEntry> matrixPositions = new LinkedHashMap();
    Map<String, Long> blockIndexPositions = new LinkedHashMap();
    Map<String, IndexEntry[]> blockIndexMap = new LinkedHashMap();


    private int countThreshold = 0;
    private boolean diagonalsOnly = false;
    private Set<String> includedChromosomes = null;

    //static DensityCalculation densityCalculation;

    public Preprocessor(File outputFile, List<Chromosome> chromosomes) {
        this.outputFile = outputFile;
        this.chromosomes = chromosomes;

        chromosomeOrdinals = new Hashtable();
        for (int i = 0; i < chromosomes.size(); i++) {
            chromosomeOrdinals.put(chromosomes.get(i).getName(), i);
        }
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

    public void preprocess(List<String> inputFileList) throws IOException {


        try {
            fos = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

            //densityCalculation = new DensityCalculation(HiCTools.chromosomes);

            // Placeholder for master index position, replace later
            writeLong(0l);

            int nChrs = chromosomes.size();
            writeInt(nChrs);
            for (Chromosome chromosome : chromosomes) {
                writeString(chromosome.getName());
                writeInt(chromosome.getSize());
            }

            // Attribute dictionary -- nothing for now, reserve for future.
            int nAttributes = 0;
            writeInt(nAttributes);
            // Future -- loop through attributes writing key/value pairs


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
                        System.out.println("writing matrix: " + matrix.getKey());
                        writeMatrix(matrix);
                    }
                }
            }


            masterIndexPosition = bytesWritten;
            writeMasterIndex();


        } finally {
            fos.close();
        }

        updateIndexPositions();
    }


    /**
     * Compute matrix for the given chromosome combination.  This resultes in full pass through the input files
     * for each chromosome combination.  This is done to save memory, at the expense of longer running times.
     *
     * @param inputFileList
     * @param c1
     * @param c2
     * @return
     * @throws IOException
     */
    public MatrixPP computeMatrix(List<String> inputFileList, int c1, int c2) throws IOException {

        boolean isWholeGenome = (c1 == 0 && c2 == 0);

        MatrixPP matrix = null;

        if (isWholeGenome) {
            int genomeLength = chromosomes.get(0).getSize();  // <= whole genome in KB
            int binSize = genomeLength / 500;
            matrix = new MatrixPP(c1, c2, binSize);
        } else {
            matrix = new MatrixPP(c1, c2);
        }

        for (String file : inputFileList) {

            PairIterator iter = (file.endsWith(".bam") || file.endsWith(".bam.hg19")) ?
                    new BAMPairIterator(file) :
                    new AsciiPairIterator(file);

            while (iter.hasNext()) {

                AlignmentPair pair = iter.next();
                int pos1 = pair.getPos1();
                int pos2 = pair.getPos2();
                Integer chr1 = chromosomeOrdinals.get(pair.getChr1());
                Integer chr2 = chromosomeOrdinals.get(pair.getChr2());
                if (chr1 != null && chr2 != null) {
                    if (isWholeGenome) {
                        pos1 = getGenomicPosition(chr1, pos1);
                        pos2 = getGenomicPosition(chr2, pos2);
                        incrementCount(matrix, c1, pos1, c2, pos2);
                        totalCount++;
                    } else if ((c1 == chr1 && c2 == chr2) || (c1 == chr2 && c2 == chr1)) {
                        incrementCount(matrix, chr1, pos1, chr2, pos2);
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
            len += chromosomes.get(i).getSize();
        }
        len += pos;

        return (int) (len / 1000);

    }

    private static void incrementCount(MatrixPP matrix, int chr1, int pos1, int chr2, int pos2) {

        if (chr2 > chr1) {
            //transpose
            int tc2 = chr2;
            int tp2 = pos2;
            chr2 = chr1;
            pos2 = pos1;
            chr1 = tc2;
            pos1 = tp2;
        }

        matrix.incrementCount(pos1, pos2);

        if (chr1 == chr2) {
            int dist = Math.abs(pos1 - pos2);
            //densityCalculation.addDistance(chr1, dist);
        }
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

            // Block indeces
            for (String key : blockIndexPositions.keySet()) {
                long pos = blockIndexPositions.get(key);
                IndexEntry[] blockIndex = blockIndexMap.get(key);

                raf.getChannel().position(pos);

                // Write as little endian
                buffer = new BufferedByteWriter();
                for (int i = 0; i < blockIndex.length; i++) {
                    buffer.putInt(blockIndex[i].id);
                    buffer.putLong(blockIndex[i].position);
                    buffer.putInt(blockIndex[i].size);
                }
                raf.write(buffer.getBytes());
            }
        } finally {
            if (raf != null) raf.close();
        }
    }

    public void writeMasterIndex() throws IOException {

        BufferedByteWriter buffer = new BufferedByteWriter();
        buffer.putInt(matrixPositions.size());
        for (Map.Entry<String, IndexEntry> entry : matrixPositions.entrySet()) {
            buffer.putString(entry.getKey());
            buffer.putLong(entry.getValue().position);
            buffer.putInt(entry.getValue().size);
        }

        //writeExpectedValues(buffer);

        byte[] bytes = buffer.getBytes();
        writeInt(bytes.length);
        write(bytes);
    }
//
//    private void writeExpectedValues(BufferedByteWriter buffer) throws IOException {
//        buffer.putLong(totalCount);
//
//        densityCalculation.computeDensity();
//
//        buffer.putInt(densityCalculation.getGridSize());
//        double[] densities = densityCalculation.getDensityAvg();
//        buffer.putInt(densities.length);
//        for(int i=0; i<densities.length; i++) {
//            buffer.putDouble(densities[i]);
//        }
//
//        Map<Integer, Double> normFactors = densityCalculation.getNormalizationFactors();
//        buffer.putInt(normFactors.size());
//        for(Map.Entry<Integer, Double> entry : normFactors.entrySet()) {
//            buffer.putInt(entry.getKey());
//            buffer.putDouble(entry.getValue());
//        }
//    }


    public void writeMatrix(MatrixPP matrix) throws IOException {

        long position = bytesWritten;

        writeInt(matrix.getChr1());
        writeInt(matrix.getChr2());
        writeInt(matrix.getZoomData().length);
        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            writeZoomHeader(zd);
        }
        int size = (int) (bytesWritten - position);
        matrixPositions.put(matrix.getKey(), new IndexEntry(position, size));


        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            IndexEntry[] blockIndex = writeZoomData(zd);
            blockIndexMap.put(getBlockKey(zd), blockIndex);
        }
    }

    private String getBlockKey(MatrixZoomDataPP zd) {
        return zd.getChr1() + "_" + zd.getChr2() + "_" + zd.getZoom();
    }

    private void writeZoomHeader(MatrixZoomDataPP zd) throws IOException {

        writeInt(zd.getZoom());
        writeInt(zd.getBinSize());
        writeInt(zd.getBlockBinCount());
        writeInt(zd.getBlockColumnCount());

        final Map<Integer, Block> blocks = zd.getBlocks();
        writeInt(blocks.size());
        blockIndexPositions.put(getBlockKey(zd), bytesWritten);

        // Placeholder for block index
        for (int i = 0; i < zd.getBlocks().size(); i++) {
            writeInt(0);
            writeLong(0l);
            writeInt(0);
        }

    }

    private IndexEntry[] writeZoomData(MatrixZoomDataPP zd) throws IOException {

        final Map<Integer, Block> blocks = zd.getBlocks();

        IndexEntry[] indexEntries = new IndexEntry[blocks.size()];
        int i = 0;
        for (Map.Entry<Integer, Block> entry : blocks.entrySet()) {

            int blockNumber = entry.getKey().intValue();
            Block block = entry.getValue();

            long position = bytesWritten;
            writeContactRecords(block);
            int size = (int) (bytesWritten - position);

            indexEntries[i] = new IndexEntry(blockNumber, position, size);
            i++;
        }
        return indexEntries;

    }

    /**
     * Note -- compressed
     *
     * @param block
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
                buffer.putInt(rec.getX());
                buffer.putInt(rec.getY());
                buffer.putInt(rec.getCounts());
            }
        }

        byte[] bytes = buffer.getBytes();
        byte[] compressedBytes = CompressionUtils.compress(bytes);
        write(compressedBytes);

    }

    private void writeDouble(double d) throws IOException {
        fos.writeDouble(d);
        bytesWritten += 8;
    }


    private void writeInt(int v) throws IOException {
        fos.writeInt(v);
        bytesWritten += 4;
    }

    private void writeShort(short v) throws IOException {

        fos.writeShort(v);
        bytesWritten += 2;
    }

    public void writeLong(long v) throws IOException {
        fos.writeLong(v);
        bytesWritten += 8;
    }

    private void write(byte[] bytes) throws IOException {
        fos.write(bytes);
        bytesWritten += bytes.length;
    }

    private void writeString(String string) throws IOException {
        byte[] bytes = string.getBytes();
        write(bytes);
        fos.write((byte) 0);
        bytesWritten++;
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
     * @date Aug 12, 2010
     */
    class MatrixPP {

        private int chr1;
        private int chr2;
        private MatrixZoomDataPP[] zoomData;


        /**
          * Constructor for creating a matrix and initializing zoomd data at predefined resolution scales.  This
          * constructor is used when parsing alignment files.
          *
          * @param chr1
          * @param chr2
          */
         MatrixPP(int chr1, int chr2) {
             this.chr1 = chr1;
             this.chr2 = chr2;
             zoomData = new MatrixZoomDataPP[HiCGlobals.zoomBinSizes.length];
             for (int zoom = 0; zoom < HiCGlobals.zoomBinSizes.length; zoom++) {
                 int binSize = HiCGlobals.zoomBinSizes[zoom];
                 int nColumns = (int) Math.pow(Math.pow(2, zoom), 0.25);
                 zoomData[zoom] = new MatrixZoomDataPP(chr1, chr2, binSize, nColumns, zoom);
             }
         }

        /**
         * Constructor for creating a matrix with a single zoom level at a specified bin size.  This is provided
         * primarily for constructing a whole-genome view.
         *
         * @param chr1
         * @param chr2
         * @param binSize
         */
        MatrixPP(int chr1, int chr2, int binSize) {
            this.chr1 = chr1;
            this.chr2 = chr2;
            zoomData = new MatrixZoomDataPP[1];
            int nBlocks = 1;
            zoomData[0] = new MatrixZoomDataPP(chr1, chr2, binSize, nBlocks, 0);

        }


        String generateKey(int chr1, int chr2) {
            return "" + chr1 + "_" + chr2;
        }

        String getKey() {
            return generateKey(chr1, chr2);
        }


        void incrementCount(int pos1, int pos2) {

            for (int i = 0; i < zoomData.length; i++) {
                zoomData[i].incrementCount(pos1, pos2);
            }

        }

        void parsingComplete() {
            for (MatrixZoomDataPP zd : zoomData) {
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
     * @date Aug 10, 2010
     */
    class MatrixZoomDataPP {

        private int chr1;  // Redundant, but convenient    BinDatasetReader
        private int chr2;  // Redundant, but convenient

        private int zoom;
        private int binSize;         // bin size in bp
        private int blockBinCount;   // block size in bins
        private int blockColumnCount;     // number of block columns

        private LinkedHashMap<Integer, Block> blocks;
        private Map<Integer, Preprocessor.IndexEntry> blockIndex;


        int getBinSize() {
            return binSize;
        }


        int getChr1() {
            return chr1;
        }


        int getChr2() {
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
         *
         * @param chr1             index of first chromosome  (x-axis)
         * @param chr2
         * @param binSize          size of each grid bin in bp
         * @param blockColumnCount number of block columns
         * @param zoom             integer zoom (resolution) level index.  TODO Is this needed?
         */
        MatrixZoomDataPP(int chr1, int chr2, int binSize, int blockColumnCount, int zoom) {


            this.chr1 = chr1;
            this.chr2 = chr2;
            this.binSize = binSize;
            this.blockColumnCount = blockColumnCount;
            this.zoom = zoom;


            int nBinsX = chromosomes.get(chr1).getSize() / binSize + 1;
            blockBinCount = nBinsX / blockColumnCount + 1;
            blocks = new LinkedHashMap(blockColumnCount * blockColumnCount);
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
