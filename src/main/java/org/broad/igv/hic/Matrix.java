package org.broad.igv.hic;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;
import org.broad.igv.ucsc.twobit.UnsignedByteBufferImpl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Matrix {

    private final int chr1;
    private final int chr2;
    private final List<MatrixZoomData> bpZoomData = new ArrayList<>();
    private final List<MatrixZoomData> fragZoomData = new ArrayList<>();

    public Matrix(int chr1, int chr2, List<MatrixZoomData> zoomDataList) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        for (MatrixZoomData zd : zoomDataList) {
            if ("BP".equals(zd.getZoom().unit())) {
                this.bpZoomData.add(zd);
            } else {
                this.fragZoomData.add(zd);
            }
        }
    }

    /**
     * Find the best zoom level for the given bin size
     */
    public int findZoomForResolution(int binSize, String unit) {
        List<MatrixZoomData> zdArray = "FRAG".equals(unit) ? this.fragZoomData : this.bpZoomData;
        for (int i = 1; i < zdArray.size(); i++) {
            MatrixZoomData zd = zdArray.get(i);
            if (zd.getZoom().binSize() < binSize) {
                return i - 1;
            }
        }
        return Math.max(0, zdArray.size() - 1);
    }

    /**
     * Fetch zoom data by bin size. Returns null if not found.
     */
    public MatrixZoomData getZoomData(int binSize, String unit) {
        if (unit == null) unit = "BP";
        List<MatrixZoomData> zdArray = "BP".equals(unit) ? this.bpZoomData : this.fragZoomData;
        for (MatrixZoomData zd : zdArray) {
            if (binSize == zd.getZoom().binSize()) {
                return zd;
            }
        }
        return null;
    }

    /**
     * Return zoom data by resolution index.
     */
    public MatrixZoomData getZoomDataByIndex(int index, String unit) {
        List<MatrixZoomData> zdArray = "FRAG".equals(unit) ? this.fragZoomData : this.bpZoomData;
        return zdArray.get(index);
    }

    public int getChr1() {
        return chr1;
    }

    public int getChr2() {
        return chr2;
    }

    public List<MatrixZoomData> getBpZoomData() {
        return bpZoomData;
    }

    public List<MatrixZoomData> getFragZoomData() {
        return fragZoomData;
    }

    public static String getKey(int chrIdx1, int chrIdx2) {
        if (chrIdx1 > chrIdx2) {
            int tmp = chrIdx1;
            chrIdx1 = chrIdx2;
            chrIdx2 = tmp;
        }
        return chrIdx1 + "_" + chrIdx2;
    }

    /**
     * Parse a matrix from binary data. Expects a UnsignedByteBuffer with methods used below.
     * Adjust parameters to match your project's binary parser (e.g. ByteBuffer, DataInputStream).
     */
    public static Matrix parseMatrix(byte[] data, List<Chromosome> chromosomes) throws IOException {
        UnsignedByteBuffer dis =  UnsignedByteBufferImpl.wrap(data); // adjust to your UnsignedByteBuffer constructor
        int c1 = dis.getInt();
        int c2 = dis.getInt();

        // optional validation could be added here
        Chromosome chr1 = chromosomes.get(c1);
        Chromosome chr2 = chromosomes.get(c2);

        int nResolutions = dis.getInt();
        List<MatrixZoomData> zdList = new ArrayList<>(nResolutions);

        while (nResolutions-- > 0) {
            MatrixZoomData zd = MatrixZoomData.parseMatrixZoomData(chr1, chr2, dis);
            zdList.add(zd);
        }
        return new Matrix(c1, c2, zdList);
    }
}