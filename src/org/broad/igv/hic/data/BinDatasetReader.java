package org.broad.igv.hic.data;


import org.broad.igv.util.CompressionUtils;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.SeekableFileStream;
import org.broad.tribble.util.SeekableStream;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Aug 17, 2010
 */
public class BinDatasetReader implements DatasetReader {

    SeekableStream stream;
    File file;

    Map<String, DatasetWriter.IndexEntry> masterIndex = new HashMap();

    public BinDatasetReader(File file) {
        this.file = file;
    }

    public Dataset read() throws FileNotFoundException {

        stream = new SeekableFileStream(file);

        try {
            long masterIndexPos = readLong();
            readMasterIndex(masterIndexPos);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        Dataset ds = new Dataset(this);

        return ds;

    }


    private Map<String, DatasetWriter.IndexEntry> readMasterIndex(long position) throws IOException {

        stream.seek(position);
        int nBytes = readInt();

        byte[] buffer = new byte[nBytes];
        stream.readFully(buffer);

        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
        int nEntries = dis.readInt();
        for (int i = 0; i < nEntries; i++) {
            String key = dis.readString();
            long filePosition = dis.readLong();
            int sizeInBytes = dis.readInt();
            masterIndex.put(key, new DatasetWriter.IndexEntry(filePosition, sizeInBytes));
        }

        return masterIndex;

    }

    public Matrix readMatrix(String key) throws IOException {
        DatasetWriter.IndexEntry idx = masterIndex.get(key);
        if (idx == null) {
            return null;
        }

        byte[] buffer = new byte[idx.size];
        stream.seek(idx.position);
        stream.readFully(buffer);

        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
        int c1 = dis.readInt();
        int c2 = dis.readInt();
        int nZooms = dis.readInt();


//    public MatrixZoomData(int chr1, int chr2, int binSize, int blockColumnCount, int zoom) {

        MatrixZoomData[] zd = new MatrixZoomData[nZooms];
        for (int i = 0; i < nZooms; i++) {
            int zoom = dis.readInt();
            int binSize = dis.readInt();
            int blockSize = dis.readInt();
            int blockColumnCount = dis.readInt();


            int nBlocks = dis.readInt();
            Map<Integer, DatasetWriter.IndexEntry> blockIndex = new HashMap(nBlocks);

            for (int b = 0; b < nBlocks; b++) {
                int blockNumber = dis.readInt();
                long filePosition = dis.readLong();
                int blockSizeInBytes = dis.readInt();
                blockIndex.put(blockNumber, new DatasetWriter.IndexEntry(filePosition, blockSizeInBytes));
           }

            zd[i] = new MatrixZoomData(c1, c2, binSize, blockSize, blockColumnCount, zoom, blockIndex, this);

        }
        return new Matrix(c1, c2, zd);
    }

    public Block readBlock(int blockNumber, DatasetWriter.IndexEntry idx) throws IOException {

        byte[] compressedBytes = new byte[idx.size];
        stream.seek(idx.position);
        stream.readFully(compressedBytes);

        byte [] buffer = CompressionUtils.decompress(compressedBytes);
        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));

        int nRecords = dis.readInt();
        ContactRecord[] records = new ContactRecord[nRecords];
        for (int i = 0; i < nRecords; i++) {
            int bin1 = dis.readInt();
            int bin2 = dis.readInt();
            short counts = dis.readShort();

            records[i] = new ContactRecord(blockNumber, bin1, bin2, counts);
        }

        return new Block(blockNumber, records);

    }


    private final int readInt() throws IOException {
        byte[] buffer = new byte[4];
        stream.readFully(buffer);
        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
        return dis.readInt();
    }

    private final long readLong() throws IOException {
        byte[] buffer = new byte[8];
        stream.readFully(buffer);
        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
        return dis.readLong();
    }


}
