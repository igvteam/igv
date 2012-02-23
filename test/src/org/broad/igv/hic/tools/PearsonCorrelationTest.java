package org.broad.igv.hic.tools;

import org.apache.commons.math.linear.OpenMapRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 2/23/12
 */
public class PearsonCorrelationTest {

    @Test
    public void test1() {
        double[][] data = {{1, 2, 3}, {2.2, 4.1, 5.5}, {-1.5, -2.5, -2.5}};

        double[][] expectedValues = {{
                1.0, 0.9961956797340459, -0.8660254037844386},
                {0.9961956797340459, 1.0, -0.9063030266884746},
                {-0.8660254037844386, -0.9063030266884746, 1.0}};    // correl(1,2)

        RealMatrix rm = new OpenMapRealMatrix(3, 3);
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                rm.addToEntry(row, col, data[col][row]);
            }
        }

        PearsonsCorrelation corr = (new PearsonsCorrelation());
        System.out.println(PearsonCorrelation.getPearsonCorrelation(data[0], data[1]));
        System.out.println(PearsonCorrelation.getPearsonCorrelation(data[0], data[2]));
        System.out.println(PearsonCorrelation.getPearsonCorrelation(data[1], data[2]));

        RealMatrix m = corr.computeCorrelationMatrix(rm);
        for (int i = 0; i < m.getRowDimension(); i++) {
            for (int j = 0; j < m.getColumnDimension(); j++) {
                double tolerance = 1.0e-6;
                assertEquals(expectedValues[i][j], m.getEntry(i, j), tolerance);
            }
            System.out.println();
        }
    }
}
