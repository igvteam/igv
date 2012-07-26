package org.broad.igv.hic.matrix;

import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.LittleEndianInputStream;

import java.io.*;

/**
 * Class for testing/development
 * <p/>
 * Assumptions -- matrix is square
 *
 * @author jrobinso
 *         Date: 7/13/12
 *         Time: 1:48 PM
 */
public class InMemoryMatrix implements BasicMatrix {


    int dim;
    float[] data;

    public InMemoryMatrix(float[] data, int dim) {
        this.data = data;
        this.dim = dim;
    }


    @Override
    public float getEntry(int row, int col) {

        int idx = row * dim + col;
        return data[idx];
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
        System.arraycopy(data, 0, subdata, 0, subdata.length);

        return new InMemoryMatrix(subdata, dim);
    }

    @Override
    public float getLowerValue() {
        return -1;
    }

    @Override
    public float getUpperValue() {
        return 1;
    }
}
