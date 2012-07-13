package org.broad.igv.hic.data;

import org.apache.commons.math.linear.*;

/**
 * Implementation of a disk-resident (as opposed to memory-resident) matrix.   This is not a full implementation,
 * only methods required for Hi-C operations are implemented.
 *
 * @author jrobinso
 *         Date: 7/13/12
 *         Time: 10:20 AM
 */
public class BlockIndexedRealMatrix implements RealMatrix{


    @Override
    public RealMatrix getSubMatrix(int i, int i1, int i2, int i3) throws MatrixIndexException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix getSubMatrix(int[] ints, int[] ints1) throws MatrixIndexException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix createMatrix(int i, int i1) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix copy() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix add(RealMatrix realMatrix) throws IllegalArgumentException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix subtract(RealMatrix realMatrix) throws IllegalArgumentException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix scalarAdd(double v) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix scalarMultiply(double v) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix multiply(RealMatrix realMatrix) throws IllegalArgumentException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix preMultiply(RealMatrix realMatrix) throws IllegalArgumentException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double[][] getData() {
        return new double[0][];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getNorm() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getFrobeniusNorm() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void copySubMatrix(int i, int i1, int i2, int i3, double[][] doubles) throws MatrixIndexException, IllegalArgumentException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void copySubMatrix(int[] ints, int[] ints1, double[][] doubles) throws MatrixIndexException, IllegalArgumentException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setSubMatrix(double[][] doubles, int i, int i1) throws MatrixIndexException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix getRowMatrix(int i) throws MatrixIndexException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setRowMatrix(int i, RealMatrix realMatrix) throws MatrixIndexException, InvalidMatrixException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix getColumnMatrix(int i) throws MatrixIndexException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setColumnMatrix(int i, RealMatrix realMatrix) throws MatrixIndexException, InvalidMatrixException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealVector getRowVector(int i) throws MatrixIndexException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setRowVector(int i, RealVector realVector) throws MatrixIndexException, InvalidMatrixException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealVector getColumnVector(int i) throws MatrixIndexException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setColumnVector(int i, RealVector realVector) throws MatrixIndexException, InvalidMatrixException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double[] getRow(int i) throws MatrixIndexException {
        return new double[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setRow(int i, double[] doubles) throws MatrixIndexException, InvalidMatrixException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double[] getColumn(int i) throws MatrixIndexException {
        return new double[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setColumn(int i, double[] doubles) throws MatrixIndexException, InvalidMatrixException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getEntry(int i, int i1) throws MatrixIndexException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setEntry(int i, int i1, double v) throws MatrixIndexException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void addToEntry(int i, int i1, double v) throws MatrixIndexException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void multiplyEntry(int i, int i1, double v) throws MatrixIndexException {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix transpose() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix inverse() throws InvalidMatrixException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getDeterminant() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean isSingular() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getTrace() throws NonSquareMatrixException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double[] operate(double[] doubles) throws IllegalArgumentException {
        return new double[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealVector operate(RealVector realVector) throws IllegalArgumentException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double[] preMultiply(double[] doubles) throws IllegalArgumentException {
        return new double[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealVector preMultiply(RealVector realVector) throws IllegalArgumentException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInRowOrder(RealMatrixChangingVisitor realMatrixChangingVisitor) throws MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInRowOrder(RealMatrixPreservingVisitor realMatrixPreservingVisitor) throws MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInRowOrder(RealMatrixChangingVisitor realMatrixChangingVisitor, int i, int i1, int i2, int i3) throws MatrixIndexException, MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInRowOrder(RealMatrixPreservingVisitor realMatrixPreservingVisitor, int i, int i1, int i2, int i3) throws MatrixIndexException, MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInColumnOrder(RealMatrixChangingVisitor realMatrixChangingVisitor) throws MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInColumnOrder(RealMatrixPreservingVisitor realMatrixPreservingVisitor) throws MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInColumnOrder(RealMatrixChangingVisitor realMatrixChangingVisitor, int i, int i1, int i2, int i3) throws MatrixIndexException, MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInColumnOrder(RealMatrixPreservingVisitor realMatrixPreservingVisitor, int i, int i1, int i2, int i3) throws MatrixIndexException, MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInOptimizedOrder(RealMatrixChangingVisitor realMatrixChangingVisitor) throws MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInOptimizedOrder(RealMatrixPreservingVisitor realMatrixPreservingVisitor) throws MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInOptimizedOrder(RealMatrixChangingVisitor realMatrixChangingVisitor, int i, int i1, int i2, int i3) throws MatrixIndexException, MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double walkInOptimizedOrder(RealMatrixPreservingVisitor realMatrixPreservingVisitor, int i, int i1, int i2, int i3) throws MatrixIndexException, MatrixVisitorException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double[] solve(double[] doubles) throws IllegalArgumentException, InvalidMatrixException {
        return new double[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RealMatrix solve(RealMatrix realMatrix) throws IllegalArgumentException, InvalidMatrixException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean isSquare() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getRowDimension() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getColumnDimension() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
