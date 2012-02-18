package org.broad.igv.hic.tools;

import org.broad.igv.hic.data.Block;
import org.broad.igv.hic.data.ContactRecord;
import org.broad.igv.hic.data.Matrix;
import org.broad.igv.hic.data.MatrixZoomData;
import org.broad.igv.util.CompressionUtils;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 * @date Aug 16, 2010
 */
public class Preprocessor {


    File outputFile;
    LittleEndianOutputStream fos;
    long bytesWritten = 0;
    long totalCount = 0;

    long masterIndexPosition;
    Map<String, IndexEntry> matrixPositions = new LinkedHashMap();
    Map<String, Long> blockIndexPositions = new LinkedHashMap();
    Map<String, IndexEntry[]> blockIndexMap = new LinkedHashMap();
    private int countThreshold = 3;

    private boolean diagonalsOnly = true;
    //static DensityCalculation densityCalculation;

    public Preprocessor(File outputFile) {
        this.outputFile = outputFile;
    }


    public void preprocess(List<String> inputFileList, String genomeId) throws IOException {


        try {
            fos = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

            //densityCalculation = new DensityCalculation(HiCTools.chromosomes);

            // Placeholder for master index position, replace later
            writeLong(0l);

            int nChrs = HiCTools.chromosomes.length;
            writeInt(nChrs);
            for (int i = 0; i < nChrs; i++) {
                writeString(HiCTools.chromosomes[i].getName());
                writeInt(HiCTools.chromosomes[i].getSize());
            }

            // Attribute dictionary -- nothing for now, reserve for future.
            int nAttributes = 0;
            writeInt(nAttributes);
            // Future -- loop through attributes writing key/value pairs


            // Data
            for (int c1 = 0; c1 < nChrs; c1++) {
                for (int c2 = c1; c2 < nChrs; c2++) {

                    // Index zero is whole genome
                    if ((c1 == 0 && c2 != 0) || (c2 == 0 && c1 != 0)) continue;

                    if(diagonalsOnly && c1 != c2) continue;

                    Matrix matrix = computeMatrix(inputFileList, c1, c2);

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
    public Matrix computeMatrix(List<String> inputFileList, int c1, int c2) throws IOException {

        boolean isWholeGenome = (c1 == 0 && c2 == 0);

        Matrix matrix = null;

        if (isWholeGenome) {
            int genomeLength = HiCTools.chromosomes[0].getSize();  // <= whole genome in KB
            int binSize = genomeLength / 500;
            matrix = new Matrix(c1, c2, binSize);
        } else {
            matrix = new Matrix(c1, c2);
        }

        for (String file : inputFileList) {

            PairIterator iter = (file.endsWith(".bam") || file.endsWith(".bam.hg19")) ?
                    new BAMPairIterator(file) :
                    new AsciiPairIterator(file);

            while (iter.hasNext()) {

                AlignmentPair pair = iter.next();
                int pos1 = pair.getPos1();
                int pos2 = pair.getPos2();
                Integer chr1 = HiCTools.chromosomeOrdinals.get(pair.getChr1());
                Integer chr2 = HiCTools.chromosomeOrdinals.get(pair.getChr2());
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

    private static int getGenomicPosition(int chr, int pos) {
        long len = 0;
        for (int i = 1; i < chr; i++) {
            len += HiCTools.chromosomes[i].getSize();
        }
        len += pos;

        return (int) (len / 1000);

    }

    private static void incrementCount(Matrix matrix, int chr1, int pos1, int chr2, int pos2) {

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


    public void writeMatrix(Matrix matrix) throws IOException {

        long position = bytesWritten;

        writeInt(matrix.chr1);
        writeInt(matrix.chr2);
        writeInt(matrix.zoomData.length);
        for (MatrixZoomData zd : matrix.zoomData) {
            writeZoomHeader(zd);
        }
        int size = (int) (bytesWritten - position);
        matrixPositions.put(matrix.getKey(), new IndexEntry(position, size));


        for (MatrixZoomData zd : matrix.zoomData) {
            IndexEntry[] blockIndex = writeZoomData(zd);
            blockIndexMap.put(getBlockKey(zd), blockIndex);
        }
    }

    private String getBlockKey(MatrixZoomData zd) {
        return zd.getChr1() + "_" + zd.getChr2() + "_" + zd.getZoom();
    }

    private void writeZoomHeader(MatrixZoomData zd) throws IOException {

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

    private IndexEntry[] writeZoomData(MatrixZoomData zd) throws IOException {

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
        final int len = records.size();

        BufferedByteWriter buffer = new BufferedByteWriter(len * 12);

        buffer.putInt(len);
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

}
