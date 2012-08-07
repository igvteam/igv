package org.broad.igv.hic.matrix;

import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.util.collections.DoubleArrayList;

/**
 * @author jrobinso
 *         Date: 7/13/12
 *         Time: 1:02 PM
 */
public class RealMatrixWrapper implements BasicMatrix {

    RealMatrix matrix;
    float lowerValue = -1;
    float upperValue = 1;

    public RealMatrixWrapper(RealMatrix matrix) {
        this.matrix = matrix;
        computePercentiles();
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

    public double[] getRow(int index) {
        return matrix.getRow(index);
    }

    public BasicMatrix getSubMatrix(int startRow, int endRow, int startCol, int endCol) {
        RealMatrix subMatrix = matrix.getSubMatrix(startRow, endRow, startCol, endCol);
        return new RealMatrixWrapper(subMatrix);
    }

    @Override
    public float getLowerValue() {
        return lowerValue;
    }

    @Override
    public float getUpperValue() {
        return upperValue;
    }


    void computePercentiles() {

        // Statistics, other attributes
        DoubleArrayList flattenedDataList = new DoubleArrayList(matrix.getColumnDimension() * matrix.getRowDimension());
        double min = 1;
        double max = -1;
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                double value = matrix.getEntry(i, j);
                if (!Double.isNaN(value) && value != 1) {
                    min = value < min ? value : min;
                    max = value > max ? value : max;
                    flattenedDataList.add(value);
                }
            }
        }

        // Stats
        double[] flattenedData = flattenedDataList.toArray();
        lowerValue = (float) StatUtils.percentile(flattenedData, 5);
        upperValue = (float) StatUtils.percentile(flattenedData, 95);
        System.out.println(lowerValue + "  " + upperValue);

    }


    public RealMatrix getMatrix() {
        return matrix;
    }
}
