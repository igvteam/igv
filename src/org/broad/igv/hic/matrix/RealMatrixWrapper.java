package org.broad.igv.hic.matrix;

import org.apache.commons.math.linear.RealMatrix;

/**
 * @author jrobinso
 *         Date: 7/13/12
 *         Time: 1:02 PM
 */
public class RealMatrixWrapper implements BasicMatrix {

    RealMatrix matrix;

    public RealMatrixWrapper(RealMatrix matrix) {
        this.matrix = matrix;
    }

    @Override
    public float getEntry(int row, int col) {
        return (float) matrix.getEntry(row, col);
    }

    @Override
    public int getRowDimension() {
        return matrix.getRowDimension();
    }

    @Override
    public int getColumnDimension() {
        return matrix.getColumnDimension();
    }

    public BasicMatrix getSubMatrix(int startRow, int endRow, int startCol, int endCol) {
        RealMatrix subMatrix = matrix.getSubMatrix(startRow, endRow, startCol, endCol);
        return new RealMatrixWrapper(subMatrix);
    }


}
