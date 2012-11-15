package org.broad.igv.hic.tools;

//import org.broad.igv.hic.MainWindow;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.collections.DownsampledDoubleArrayList;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.zip.Deflater;

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

    private List<Chromosome> chromosomes;

    // Map of name -> index
    private Map<String, Integer> chromosomeIndexes;

    private File outputFile;
    private LittleEndianOutputStream los;

    private long masterIndexPosition;
    private Map<String, IndexEntry> matrixPositions;

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
    private final Deflater compressor;

    File tmpDir;

    public Preprocessor(File outputFile, String genomeId, List<Chromosome> chromosomes) {
        this.genomeId = genomeId;
        this.outputFile = outputFile;
        this.matrixPositions = new LinkedHashMap<String, IndexEntry>();

        this.chromosomes = chromosomes;
        chromosomeIndexes = new Hashtable<String, Integer>();
        for (int i = 0; i < chromosomes.size(); i++) {
            chromosomeIndexes.put(chromosomes.get(i).getName(), i);
        }

        compressor = new Deflater();
        compressor.setLevel(Deflater.DEFAULT_COMPRESSION);

        this.tmpDir = null;  // TODO -- specify this

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
                fragmentCalculation = FragmentCalculation.readFragments(fragmentFileName);
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

            los = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

            writeHeader();

            writeBody(inputFileList);

            writeFooter();


        } finally {
            if (los != null)
                los.close();
        }

        updateMasterIndex();
    }

    private void writeHeader() throws IOException {
        // Magic number
        byte[] magicBytes = "HIC".getBytes();
        los.write(magicBytes[0]);
        los.write(magicBytes[1]);
        los.write(magicBytes[2]);
        los.write(0);

        // Version
        los.writeInt(4);

        // Placeholder for master index position, replaced with actual position after all contents are written
        masterIndexPositionPosition = los.getWrittenCount();
        los.writeLong(0l);


        // Genome ID
        los.writeString(genomeId);

        // Sequence dictionary
        int nChrs = chromosomes.size();
        los.writeInt(nChrs);
        for (Chromosome chromosome : chromosomes) {
            los.writeString(chromosome.getName());
            los.writeInt(chromosome.getLength());
        }

        //BP resolution levels
        int nBpRes = bpBinSizes.length;
        los.writeInt(nBpRes);
        for (int i = 0; i < nBpRes; i++) {
            los.writeInt(bpBinSizes[i]);
        }

        //fragment resolutions
        int nFragRes = fragmentCalculation == null ? 0 : fragBinSizes.length;
        los.writeInt(nFragRes);
        for (int i = 0; i < nFragRes; i++) {
            los.writeInt(fragBinSizes[i]);
        }

        // fragment sites
        if (nFragRes > 0) {
            for (Chromosome chromosome : chromosomes) {
                int[] sites = fragmentCalculation.getSites(chromosome.getName());
                int nSites = sites == null ? 0 : sites.length;
                los.writeInt(nSites);
                for (int i = 0; i < nSites; i++) {
                    los.writeInt(sites[i]);
                }
            }
        }

        //hemi-frag levels (reserved for future use)
        int nHemiFrags = 0;
        los.writeInt(nHemiFrags);


        // Attribute dictionary
        int nAttributes = 0;
        los.writeInt(nAttributes);
        //fos.writeString("theKey");
        //fos.writeString("theValue");
        //... repeat for each attribute
    }

    private void writeBody(List<String> inputFileList) throws IOException {
        int nChrs = chromosomes.size();
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
                System.out.println("Available memory: " + RuntimeUtils.getAvailableMemory());
            }
        } // End of double loop through chromosomes


        masterIndexPosition = los.getWrittenCount();
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

            PairIterator iter = null;

            try {
                iter = (file.endsWith(".bin")) ?
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
            } finally {
                if (iter != null) iter.close();
            }
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

    public void updateMasterIndex() throws IOException {
        RandomAccessFile raf = null;
        try {
            raf = new RandomAccessFile(outputFile, "rw");

            // Master index
            raf.getChannel().position(masterIndexPositionPosition);
            BufferedByteWriter buffer = new BufferedByteWriter();
            buffer.putLong(masterIndexPosition);
            raf.write(buffer.getBytes());

        } finally {
            if (raf != null) raf.close();
        }
    }

    public void writeFooter() throws IOException {

        // Index
        BufferedByteWriter buffer = new BufferedByteWriter();
        buffer.putInt(matrixPositions.size());
        for (Map.Entry<String, IndexEntry> entry : matrixPositions.entrySet()) {
            buffer.putNullTerminatedString(entry.getKey());
            buffer.putLong(entry.getValue().position);
            buffer.putInt(entry.getValue().size);
        }

        // Vectors  (Expected values,  other).
        buffer.putInt(expectedValueCalculations.size());
        for (Map.Entry<String, ExpectedValueCalculation> entry : expectedValueCalculations.entrySet()) {
            String key = entry.getKey();
            ExpectedValueCalculation ev = entry.getValue();

            ev.computeDensity();

            int binSize = ev.getGridSize();
            String unit = ev.isFrag ? "FRAG" : "BP";

            buffer.putNullTerminatedString(unit);
            buffer.putInt(binSize);

            // The density values
            double[] expectedValues = ev.getDensityAvg();
            buffer.putInt(expectedValues.length);
            for (int i = 0; i < expectedValues.length; i++) {
                buffer.putDouble(expectedValues[i]);
            }

            // Map of chromosome index -> normalization factor
            Map<Integer, Double> normalizationFactors = ev.getNormalizationFactors();
            buffer.putInt(normalizationFactors.size());
            for (Map.Entry<Integer, Double> normFactor : normalizationFactors.entrySet()) {
                buffer.putInt(normFactor.getKey());
                buffer.putDouble(normFactor.getValue());
                //System.out.println(normFactor.getKey() + "  " + normFactor.getValue());
            }

        }

        byte[] bytes = buffer.getBytes();
        los.writeInt(bytes.length);
        los.write(bytes);
    }

    public synchronized void writeMatrix(MatrixPP matrix) throws IOException {

        System.out.println("Start writing matrix: " + matrix.getKey());

        long position = los.getWrittenCount();

        los.writeInt(matrix.getChr1Idx());
        los.writeInt(matrix.getChr2Idx());
        int numResolutions = 0;

        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            if (zd != null) {
                numResolutions++;
            }
        }
        los.writeInt(numResolutions);

        //fos.writeInt(matrix.getZoomData().length);
        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            if (zd != null)
                writeZoomHeader(zd);
        }

        int size = (int) (los.getWrittenCount() - position);
        matrixPositions.put(matrix.getKey(), new IndexEntry(position, size));

        for (MatrixZoomDataPP zd : matrix.getZoomData()) {
            if (zd != null) {
                List<IndexEntry> blockIndex = zd.mergeAndWriteBlocks();
                zd.updateIndexPositions(blockIndex);
            }
        }

        System.out.println("Done writing matrix: " + matrix.getKey());
    }

    private String getBlockKey(MatrixZoomDataPP zd) {
        return zd.getChr1() + "_" + zd.getChr2() + "_" + zd.getUnit() + "_" + zd.getBinSize();
    }

    private void writeZoomHeader(MatrixZoomDataPP zd) throws IOException {

        int numberOfBlocks = zd.blockNumbers.size();

//        System.out.println("Write zoom header");
//        System.out.println(zd.getUnit());  // Unit, ether "BP" or "FRAG"
//        System.out.println(zd.getZoom());     // zoom index,  lowest res is zero
//        System.out.println((float) zd.getSum());      // sum
//        System.out.println((float) zd.getAverage());
//        System.out.println((float) zd.getStdDev());
//        System.out.println((float) zd.getPercent95());
//        System.out.println(zd.getBinSize());
//        System.out.println(zd.getBlockBinCount());
//        System.out.println(zd.getBlockColumnCount());
//        System.out.println(numberOfBlocks);

        los.writeString(zd.getUnit());  // Unit, ether "BP" or "FRAG"
        los.writeInt(zd.getZoom());     // zoom index,  lowest res is zero
        los.writeFloat((float) zd.getSum());      // sum
        los.writeFloat((float) zd.getOccupiedCellCount());
        los.writeFloat((float) zd.getPercent5());
        los.writeFloat((float) zd.getPercent95());
        los.writeInt(zd.getBinSize());
        los.writeInt(zd.getBlockBinCount());
        los.writeInt(zd.getBlockColumnCount());
        los.writeInt(numberOfBlocks);

        zd.blockIndexPosition = los.getWrittenCount();

        // Placeholder for block index
        for (int i = 0; i < numberOfBlocks; i++) {
            los.writeInt(0);
            los.writeLong(0l);
            los.writeInt(0);
        }

    }

    /**
     * Note -- compressed
     *
     * @param zd
     * @param block       Block to write
     * @param sampledData Array to hold a sample of the data (to compute statistics)
     * @throws IOException
     */
    private void writeContactRecords(MatrixZoomDataPP zd, BlockPP block, DownsampledDoubleArrayList sampledData) throws IOException {

        final Map<Point, ContactCount> records = block.getContractRecordMap();//   getContactRecords();

        // System.out.println("Write contact records : records count = " + records.size());

        // Count records first
        int nRecords;
        if (countThreshold > 0) {
            nRecords = 0;
            for (ContactCount rec : records.values()) {
                if (rec.getCounts() >= countThreshold) {
                    nRecords++;
                }
            }
        } else {
            nRecords = records.size();
        }

        zd.cellCount += nRecords;

        BufferedByteWriter buffer = new BufferedByteWriter(nRecords * 12);
        buffer.putInt(nRecords);
        for (Map.Entry<Point, ContactCount> entry : records.entrySet()) {
            Point point = entry.getKey();
            float counts = entry.getValue().getCounts();
            if (counts >= countThreshold) {
                buffer.putInt(point.x);
                buffer.putInt(point.y);
                buffer.putFloat(counts);

                sampledData.add(counts);
                zd.sum += counts;
            }
        }

        byte[] bytes = buffer.getBytes();
        byte[] compressedBytes = compress(bytes);
        los.write(compressedBytes);

    }

    public void setTmpdir(String tmpDirName) {
        this.tmpDir = new File(tmpDirName);

        if(!tmpDir.exists()) {
            System.err.println("Tmp directory does not exist: " + tmpDirName);
            System.exit(1);
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
         * @param chr1Idx Chromosome 1
         * @param chr2Idx Chromosome 2
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


        void incrementCount(int pos1, int pos2, float score) throws IOException {

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

        private double sum = 0;
        private double cellCount = 0;
        private double percent5;
        private double percent95;

        private int zoom;
        private int binSize;         // bin size in bp
        private int blockBinCount;   // block size in bins
        private int blockColumnCount;     // number of block columns

        boolean isFrag;

        private LinkedHashMap<Integer, BlockPP> blocks;


        Set<Integer> blockNumbers;  // The only reason for this is to get a count

        List<File> tmpFiles;
        public long blockIndexPosition;

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

            this.tmpFiles = new ArrayList<File>();
            this.blockNumbers = new HashSet<Integer>(1000);

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
            blocks = new LinkedHashMap<Integer, BlockPP>(blockColumnCount * blockColumnCount);
        }

        String getUnit() {
            return isFrag ? "FRAG" : "BP";
        }

        double getSum() {
            return sum;
        }

        public double getOccupiedCellCount() {
            return cellCount;
        }

        public double getPercent95() {
            return percent95;
        }

        public double getPercent5() {
            return percent5;
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

        Map<Integer, BlockPP> getBlocks() {
            return blocks;
        }

        /**
         * Increment the count for the bin represented by the GENOMIC position (pos1, pos2)
         *
         * @param gPos1
         * @param gPos2
         */
        public void incrementCount(int gPos1, int gPos2, float score) throws IOException {

            sum += score;
            // Convert to proper units,  fragments or base-pairs

            int pos1 = isFrag ? fragmentCalculation.getBin(chr1.getName(), gPos1) : gPos1;
            int pos2 = isFrag ? fragmentCalculation.getBin(chr2.getName(), gPos2) : gPos2;
            if (pos1 < 0 || pos2 < 0) return;

            int xBin = pos1 / binSize;
            int yBin = pos2 / binSize;

            // Intra chromosome -- we'll store lower diagonal only
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

            BlockPP block = blocks.get(blockNumber);
            if (block == null) {

                block = new BlockPP(blockNumber);
                blocks.put(blockNumber, block);
            }
            block.incrementCount(xBin, yBin, score);

            // If too many blocks write to tmp directory
            if (blocks.size() > 10000) {
                File tmpfile = tmpDir == null ? File.createTempFile("blocks", "bin") : File.createTempFile("blocks", "bin", tmpDir);
                tmpfile.deleteOnExit();

                System.out.println(chr1.getName() + "-" + chr2.getName() + " Dumping blocks to " + tmpfile.getAbsolutePath());

                dumpBlocks(tmpfile);
                blocks.clear();
            }
        }


        /**
         * Dump the blocks calculated so far to a temporary file
         *
         * @param file
         * @throws IOException
         */
        private void dumpBlocks(File file) throws IOException {
            LittleEndianOutputStream los = null;
            try {
                los = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(file)));

                List<BlockPP> blockList = new ArrayList<BlockPP>(blocks.values());
                Collections.sort(blockList, new Comparator<BlockPP>() {
                    @Override
                    public int compare(BlockPP o1, BlockPP o2) {
                        return o1.getNumber() - o2.getNumber();
                    }
                });

                for (int i = 0; i < blockList.size(); i++) {
                    BlockPP b = blockList.get(i);
                    int number = b.getNumber();
                    blockNumbers.add(number);

                    los.writeInt(number);
                    Map<Point, ContactCount> records = b.getContractRecordMap();

                    los.writeInt(records.size());
                    for (Map.Entry<Point, ContactCount> entry : records.entrySet()) {

                        Point point = entry.getKey();
                        ContactCount count = entry.getValue();

                        los.writeInt(point.x);
                        los.writeInt(point.y);
                        los.writeFloat(count.getCounts());
                    }
                }

            } finally {
                if (los == null) los.close();
            }
        }


        // Merge and write out blocks one at a time.
        private List<IndexEntry> mergeAndWriteBlocks() throws IOException {

            DownsampledDoubleArrayList sampledData = new DownsampledDoubleArrayList(10000);

            List<BlockQueue> activeList = new ArrayList<BlockQueue>();

            // Initialize queues -- first whatever is left over in memory
            if (blocks.size() > 0) {
                BlockQueue bqInMem = new BlockQueueMem(blocks.values());
                activeList.add(bqInMem);
            }

            for (File file : tmpFiles) {
                BlockQueue bq = new BlockQueueFB(file);
                if (bq.getBlock() != null) activeList.add(bq);
            }

            List<IndexEntry> indexEntries = new ArrayList<IndexEntry>();

            do {
                Collections.sort(activeList, new Comparator<BlockQueue>() {
                    @Override
                    public int compare(BlockQueue o1, BlockQueue o2) {
                        return o1.getBlock().getNumber() - o2.getBlock().getNumber();
                    }
                });

                BlockQueue topQueue = activeList.get(0);
                BlockPP currentBlock = topQueue.getBlock();
                topQueue.advance();
                int num = currentBlock.getNumber();


                for (int i = 1; i < activeList.size(); i++) {
                    BlockQueue blockQueue = activeList.get(i);
                    if (blockQueue.getBlock().getNumber() == num) {
                        currentBlock.merge(blockQueue.getBlock());
                        blockQueue.advance();
                    }
                }

                Iterator<BlockQueue> iterator = activeList.iterator();
                while (iterator.hasNext()) {
                    if (iterator.next().getBlock() == null) {
                        iterator.remove();
                    }
                }

                // Output block
                long position = los.getWrittenCount();
                writeContactRecords(this, currentBlock, sampledData);
                int size = (int) (los.getWrittenCount() - position);

                indexEntries.add(new IndexEntry(num, position, size));


            } while (activeList.size() > 0);

            for (File f : tmpFiles) {
                f.delete();
            }

            computeStats(sampledData);

            return indexEntries;
        }

        private void computeStats(DownsampledDoubleArrayList sampledData) {

            double[] data = sampledData.toArray();
            this.percent5 = StatUtils.percentile(data, 5);
            this.percent95 = StatUtils.percentile(data, 95);

        }

        public void parsingComplete() {
            // Add the block numbers still in memory
            for (BlockPP block : blocks.values()) {
                blockNumbers.add(block.getNumber());
            }
        }

        public void updateIndexPositions(List<IndexEntry> blockIndex) throws IOException {

            // Temporarily close output stream.  Remember position
            long losPos = los.getWrittenCount();
            los.close();

            RandomAccessFile raf = null;
            try {
                raf = new RandomAccessFile(outputFile, "rw");

                // Master index
                raf.getChannel().position(masterIndexPositionPosition);
                BufferedByteWriter buffer = new BufferedByteWriter();
                buffer.putLong(masterIndexPosition);
                raf.write(buffer.getBytes());

                // Block indices
                long pos = blockIndexPosition;
                raf.getChannel().position(pos);

                // Write as little endian
                buffer = new BufferedByteWriter();
                for (IndexEntry aBlockIndex : blockIndex) {
                    buffer.putInt(aBlockIndex.id);
                    buffer.putLong(aBlockIndex.position);
                    buffer.putInt(aBlockIndex.size);
                }
                raf.write(buffer.getBytes());

            } finally {

                if (raf != null) raf.close();

                // Restore
                FileOutputStream fos = new FileOutputStream(outputFile);
                fos.getChannel().position(losPos);
                los = new LittleEndianOutputStream(new BufferedOutputStream(fos));
                los.setWrittenCount(losPos);

            }
        }
    }

    private synchronized byte[] compress(byte[] data) {

        // Give the compressor the data to compress
        compressor.reset();
        compressor.setInput(data);
        compressor.finish();

        // Create an expandable byte array to hold the compressed data.
        // You cannot use an array that's the same size as the orginal because
        // there is no guarantee that the compressed data will be smaller than
        // the uncompressed data.
        ByteArrayOutputStream bos = new ByteArrayOutputStream(data.length);

        // Compress the data
        byte[] buf = new byte[1024];
        while (!compressor.finished()) {
            int count = compressor.deflate(buf);
            bos.write(buf, 0, count);
        }
        try {
            bos.close();
        } catch (IOException e) {
            System.err.println("Error clossing ByteArrayOutputStream");
            e.printStackTrace();
        }

        byte[] compressedData = bos.toByteArray();
        return compressedData;
    }


    // class to support block merging

    static interface BlockQueue {

        void advance() throws IOException;

        BlockPP getBlock();

    }

    static class BlockQueueFB implements BlockQueue {

        File file;
        BlockPP block;
        long filePosition;

        BlockQueueFB(File file) {
            this.file = file;
            try {
                advance();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        public void advance() throws IOException {

            if (filePosition >= file.length()) {
                block = null;
                return;
            }

            FileInputStream fis = new FileInputStream(file);
            fis.getChannel().position(filePosition);


            LittleEndianInputStream lis = new LittleEndianInputStream(new BufferedInputStream(fis));


            int blockNumber = lis.readInt();
            int nRecords = lis.readInt();

            Map<Point, ContactCount> contactRecordMap = new HashMap<Point, ContactCount>(nRecords);
            for (int i = 0; i < nRecords; i++) {
                int x = lis.readInt();
                int y = lis.readInt();
                float v = lis.readFloat();
                ContactCount rec = new ContactCount(v);
                contactRecordMap.put(new Point(x, y), rec);
            }
            block = new BlockPP(blockNumber, contactRecordMap);

            // Update file position based on # of bytes read, for next block
            int nBytes = (1 + 1 + nRecords * 3) * 4;   // ints and floats are 4 bytes
            filePosition += nBytes;
        }

        public BlockPP getBlock() {
            return block;
        }
    }

    static class BlockQueueMem implements BlockQueue {

        List<BlockPP> blocks;
        int idx = 0;

        BlockQueueMem(Collection<BlockPP> blockCollection) {

            this.blocks = new ArrayList<BlockPP>(blockCollection);
            Collections.sort(blocks, new Comparator<BlockPP>() {
                @Override
                public int compare(BlockPP o1, BlockPP o2) {
                    return o1.getNumber() - o2.getNumber();
                }
            });
        }

        public void advance() {
            idx++;
        }

        public BlockPP getBlock() {
            if (idx >= blocks.size()) {
                return null;
            } else {
                return blocks.get(idx);
            }
        }
    }


    /**
     * Representation of a sparse matrix block used for preprocessing.
     */
    static class BlockPP {

        private int number;

        private Map<Point, ContactCount> contactRecordMap;


        BlockPP(int number) {
            this.number = number;
            this.contactRecordMap = new HashMap<Point, ContactCount>();
        }

        public BlockPP(int number, Map<Point, ContactCount> contactRecordMap) {
            this.number = number;
            this.contactRecordMap = contactRecordMap;
        }


        public int getNumber() {
            return number;
        }

        public void incrementCount(int col, int row, float score) {
            Point p = new Point(col, row);
            ContactCount rec = contactRecordMap.get(p);
            if (rec == null) {
                rec = new ContactCount(1);
                contactRecordMap.put(p, rec);

            } else {
                rec.incrementCount(score);
            }
        }

        public void parsingComplete() {

        }

        public Map<Point, ContactCount> getContractRecordMap() {
            return contactRecordMap;
        }

        public void merge(BlockPP other) {

            for (Map.Entry<Point, ContactCount> entry : other.getContractRecordMap().entrySet()) {

                Point point = entry.getKey();
                ContactCount otherValue = entry.getValue();

                ContactCount value = contactRecordMap.get(point);
                if (value == null) {
                    contactRecordMap.put(point, otherValue);
                } else {
                    value.incrementCount(otherValue.getCounts());
                }

            }
        }
    }


    static class ContactCount {
        float value;

        ContactCount(float value) {
            this.value = value;
        }

        void incrementCount(float increment) {
            value += increment;
        }

        public float getCounts() {
            return value;
        }
    }


}
