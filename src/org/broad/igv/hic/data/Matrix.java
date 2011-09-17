package org.broad.igv.hic.data;

import org.broad.igv.hic.MainWindow;

/**
 * @author jrobinso
 * @date Aug 12, 2010
 */
public class Matrix {

    public int chr1;
    public int chr2;
    public MatrixZoomData[] zoomData;

    public Matrix(int chr1, int chr2) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        zoomData = new MatrixZoomData[MainWindow.zoomBinSizes.length];
        for (int i = 0; i < MainWindow.zoomBinSizes.length; i++) {
            int binSize = MainWindow.zoomBinSizes[i];
            int nBlocks = (int) Math.pow(Math.pow(2, i), 0.25);
            zoomData[i] = new MatrixZoomData(chr1, chr2, binSize, nBlocks, i);
        }
    }

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

    public MatrixZoomData getZoomData(int zoomIndex) {
        return zoomData[zoomIndex];
    }


    public void incrementCount(int pos1, int pos2) {

        for (int i = 0; i < zoomData.length; i++) {
            zoomData[i].incrementCount(pos1, pos2);
        }

    }


    //static int[] zoomBinSizes = {1000};

    public void parsingComplete() {
        for (MatrixZoomData zd : zoomData) {
            zd.parsingComplete();
        }
    }
}
