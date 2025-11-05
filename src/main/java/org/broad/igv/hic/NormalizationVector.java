package org.broad.igv.hic;

import htsjdk.samtools.seekablestream.SeekableStream;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;

/**
 * Converted from JavaScript NormalizationVector.
 * Reads floats/doubles from a file at a given position, caches a windowed range,
 * and returns a requested slice as a double[].
 */
public class NormalizationVector {

    public static final int DOUBLE = 8;
    public static final int FLOAT = 4;
    public static final int LONG = 8;
    public static final int INT = 4;

    private final SeekableStream fileChannel;
    private final long filePosition;
    private final int nValues;
    private final int dataType;

    // optional metadata (may be null / 0 if not provided)
    private final String type;
    private final int chrIdx;
    private final String unit;
    private final int resolution;

    private Cache cache;

    private static class Cache {
        final int start;
        final int end;
        final double[] values;

        Cache(int start, int end, double[] values) {
            this.start = start;
            this.end = end;
            this.values = values;
        }
    }

    public NormalizationVector(SeekableStream fileChannel, long filePosition, int nValues, int dataType) {
        this(fileChannel, filePosition, nValues, dataType, null, 0, null, 0);
    }

    public NormalizationVector(SeekableStream fileChannel, long filePosition, int nValues, int dataType,
                               String type, int chrIdx, String unit, int resolution) {
        this.fileChannel = fileChannel;
        this.filePosition = filePosition;
        this.nValues = nValues;
        this.dataType = dataType;
        this.type = type;
        this.chrIdx = chrIdx;
        this.unit = unit;
        this.resolution = resolution;
        this.cache = null;
    }

    /**
     * Return values in range [start, end) as a double[].
     * Returns null if the underlying read fails.
     */
    public double[] getValues(int start, int end) throws IOException {
        if (start < 0) start = 0;
        if (end > nValues) end = nValues;
        if (start >= end) return new double[0];

        if (cache == null || start < cache.start || end > cache.end) {
            int adjustedStart = Math.max(0, start - 1000);
            int adjustedEnd = Math.min(nValues, end + 1000);
            int n = adjustedEnd - adjustedStart;
            long startPosition = filePosition + (long) adjustedStart * dataType;
            int sizeInBytes = n * dataType;

            byte [] byteArray = new byte[sizeInBytes];
            int read = 0;
            while (read < sizeInBytes) {
                fileChannel.seek(startPosition + read);
                int r = fileChannel.read(byteArray, read, sizeInBytes);
                if (r < 0) break;
                read += r;
            }
            if (read < sizeInBytes) {
                return null;
            }
            ByteBuffer buf = ByteBuffer.wrap(byteArray);

            double[] values = new double[n];
            for (int i = 0; i < n; i++) {
                if (dataType == DOUBLE) {
                    values[i] = buf.getDouble();
                } else {
                    values[i] = buf.getFloat();
                }
            }
            this.cache = new Cache(adjustedStart, adjustedEnd, values);
        }

        int sliceStart = start - cache.start;
        int sliceLength = end - start;
        return Arrays.copyOfRange(cache.values, sliceStart, sliceStart + sliceLength);
    }

    public String getKey() {
        return getNormalizationVectorKey(type, chrIdx, unit, resolution);
    }

    public static String getNormalizationVectorKey(String type, int chrIdx, String unit, int resolution) {
        return type + "_" + chrIdx + "_" + unit + "_" + resolution;
    }

    // getters for fields if needed
    public long getFilePosition() {
        return filePosition;
    }

    public int getNValues() {
        return nValues;
    }

    public int getDataType() {
        return dataType;
    }
}