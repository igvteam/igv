package org.broad.igv.hic.data;

import org.broad.igv.hic.MainWindow;

/**
 * @author jrobinso
 * @date Aug 12, 2010
 */
public class Matrix {

    private int chr1;
    private int chr2;
    private MatrixZoomData[] zoomData;


    /**
     * Constructor for creating a matrix from precomputed data.
     *
     * @param chr1
     * @param chr2
     * @param zoomData
     */
    public Matrix(int chr1, int chr2, MatrixZoomData[] zoomData) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        this.zoomData = zoomData;
    }

    public static String generateKey(int chr1, int chr2) {
        return "" + chr1 + "_" + chr2;
    }

    public String getKey() {
        return generateKey(chr1, chr2);
    }

    public MatrixZoomData getObservedMatrix(int zoomIndex) {
        return zoomData[zoomIndex];
    }


}
