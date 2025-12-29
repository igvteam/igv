package org.igv.hic;

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


    private final SeekableStream fileChannel;
    private final long filePosition;
    private final int nValues;
    private final HicFile.DataType dataType;
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

    public NormalizationVector(SeekableStream fileChannel, long filePosition, int nValues, HicFile.DataType dataType) {
        this(fileChannel, filePosition, nValues, dataType, null, 0, null, 0);
    }

    public NormalizationVector(SeekableStream fileChannel, long filePosition, int nValues, HicFile.DataType dataType,
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
    public double[] getValues(int startBin, int endBin) throws IOException {

        if (startBin < 0) startBin = 0;
        if (endBin > nValues) endBin = nValues;
        if (startBin >= endBin) return new double[0];

        if (cache == null || startBin < cache.start || endBin > cache.end) {
            int adjustedStart = Math.max(0, startBin - 1000);
            int adjustedEnd = Math.min(nValues, endBin + 1000);
            int n = adjustedEnd - adjustedStart;
            long startPosition = filePosition + (long) adjustedStart * dataType.getByteSize();
            int sizeInBytes = n * dataType.getByteSize();

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
                if (dataType == HicFile.DataType.DOUBLE) {
                    values[i] = buf.getDouble();
                } else {
                    values[i] = buf.getFloat();
                }
            }
            this.cache = new Cache(adjustedStart, adjustedEnd, values);
        }

        int sliceStart = startBin - cache.start;
        int sliceLength = endBin - startBin;
        return Arrays.copyOfRange(cache.values, sliceStart, sliceStart + sliceLength);
    }

    public String getKey() {
        return getNormalizationVectorKey(type, chrIdx, unit, resolution);
    }

    public static String getNormalizationVectorKey(String type, int chrIdx, String unit, int resolution) {
        return type + "_" + chrIdx + "_" + unit + "_" + resolution;
    }

}