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

  /**
     * Constructor for creating a matrix with a single zoom level at a specified bin size.  This is provided
     * primarily for constructing a whole-genome view.
     *
     * @param chr1
     * @param chr2
     * @param binSize
     */
    public Matrix(int chr1, int chr2, int binSize) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        zoomData = new MatrixZoomData[1];
        int nBlocks = 1;
        zoomData[0] = new MatrixZoomData(chr1, chr2, binSize, nBlocks, 0);

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


     ////////////////////////////////////////////////////////////////
     // Methods below used during preprocessing only
     ////////////////////////////////////////////////////////////////

    /**
     * Constructor for creating a matrix and initializing zoomd data at predefined resolution scales.  This
     * constructor is used when parsing alignment files.
     *
     * @param chr1
     * @param chr2
     */
    public Matrix(int chr1, int chr2) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        zoomData = new MatrixZoomData[MainWindow.zoomBinSizes.length];
        for (int zoom = 0; zoom < MainWindow.zoomBinSizes.length; zoom++) {
            int binSize = MainWindow.zoomBinSizes[zoom];
            int nColumns = (int) Math.pow(Math.pow(2, zoom), 0.25);
            zoomData[zoom] = new MatrixZoomData(chr1, chr2, binSize, nColumns, zoom);
        }
    }


    public void incrementCount(int pos1, int pos2) {

        for (int i = 0; i < zoomData.length; i++) {
            zoomData[i].incrementCount(pos1, pos2);
        }

    }

    public void parsingComplete() {
        for (MatrixZoomData zd : zoomData) {
            zd.parsingComplete();
        }
    }

    public int getChr1() {
        return chr1;
    }

    public void setChr1(int chr1) {
        this.chr1 = chr1;
    }

    public int getChr2() {
        return chr2;
    }

    public void setChr2(int chr2) {
        this.chr2 = chr2;
    }

    public MatrixZoomData[] getZoomData() {
        return zoomData;
    }

    public void setZoomData(MatrixZoomData[] zoomData) {
        this.zoomData = zoomData;
    }
}
