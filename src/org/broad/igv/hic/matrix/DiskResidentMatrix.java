package org.broad.igv.hic.matrix;

import org.broad.igv.hic.HiC;
import org.broad.igv.hic.MainWindow;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ObjectCache;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.io.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Class for testing/development
 * <p/>
 * Assumptions -- matrix is square
 *
 * @author jrobinso
 *         Date: 7/13/12
 *         Time: 1:48 PM
 */
public class DiskResidentMatrix implements BasicMatrix {


    String path;
    int dim;
    float lowerValue;
    float upperValue;
    int arrayStartPosition;
    boolean isLoading = false;

    ObjectCache<Integer, float[]> rowDataCache = new ObjectCache<Integer, float[]>(10000);

    public DiskResidentMatrix(String path, int dim) throws IOException {
        this(path, 4, dim, -1, 1);
    }

    public DiskResidentMatrix(String path, int arrayStartPosition, int dim, float lowerValue, float upperValue) throws IOException {
        this.dim = dim;
        this.path = path;
        this.lowerValue = lowerValue;
        this.upperValue = upperValue;
        this.arrayStartPosition = arrayStartPosition;
    }

    public float getLowerValue() {
        return lowerValue;
    }

    public float getUpperValue() {
        return upperValue;
    }

    public void loadRowData(int startRow, int lastRow) {
        //System.out.println("Loading row " + row);
        SeekableStream is = null;
        try {
            is = SeekableStreamFactory.getStreamFor(path);

            int startFilePosition = arrayStartPosition + (startRow * dim) * 4;
            int nBytes = (lastRow - startRow) * dim * 4;
            byte[] byteArray = new byte[nBytes];

            is.seek(startFilePosition);
            is.readFully(byteArray);

            ByteArrayInputStream bis = new ByteArrayInputStream(byteArray);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);


            for (int row = startRow; row < lastRow; row++) {
                float[] rowData = new float[dim];
                for (int i = 0; i < dim; i++) {
                    float f = les.readFloat();
                    rowData[i] = f;
                }
                rowDataCache.put(row, rowData);
            }

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (is != null)
                try {
                    is.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
        }
    }

    @Override
    public float getEntry(final int row, int col) {

        if (!rowDataCache.containsKey(row)) {
            MainWindow.getInstance().showGlassPane();
            int lastRow = Math.min(dim, row + 500);
            loadRowData(row, lastRow);
            MainWindow.getInstance().hideGlassPane();
        }
        float[] rowData = rowDataCache.get(row);
        return rowData == null ? Float.NaN : rowData[col];
    }

    @Override
    public int getRowDimension() {
        return dim;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getColumnDimension() {
        return dim;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public BasicMatrix getSubMatrix(int startRow, int endRow, int startCol, int endCol) {

        int startIdx = startRow * dim + startCol;
        int endIdx = endRow * dim + endCol;
        int dim = endRow - startRow + 1;
        // TODO -- assert matrix is square
        float[] subdata = new float[endIdx - startIdx + 1];

        for (int row = startRow; row <= endRow; row++) {
            for (int col = startCol; col <= endCol; col++) {
                int idx = (row - startRow) * dim + row + (col - startCol) + col;
                subdata[idx] = getEntry(row, col);
            }
        }

        return new InMemoryMatrix(subdata, dim);
    }

}
