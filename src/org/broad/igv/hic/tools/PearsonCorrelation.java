package org.broad.igv.hic.tools;

import org.apache.commons.math.linear.OpenMapRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Jim Robinson
 * @date 12/19/11
 */
public class PearsonCorrelation {


    public static double getPearsonCorrelation(double[] scores1, double[] scores2) {
        double result = 0;
        double sum_sq_x = 0;
        double sum_sq_y = 0;
        double sum_coproduct = 0;
        double mean_x = scores1[0];
        double mean_y = scores2[0];
        for (int i = 2; i < scores1.length + 1; i += 1) {
            double sweep = Double.valueOf(i - 1) / i;
            double delta_x = scores1[i - 1] - mean_x;
            double delta_y = scores2[i - 1] - mean_y;
            sum_sq_x += delta_x * delta_x * sweep;
            sum_sq_y += delta_y * delta_y * sweep;
            sum_coproduct += delta_x * delta_y * sweep;
            mean_x += delta_x / i;
            mean_y += delta_y / i;
        }
        double pop_sd_x = (double) Math.sqrt(sum_sq_x / scores1.length);
        double pop_sd_y = (double) Math.sqrt(sum_sq_y / scores1.length);
        double cov_x_y = sum_coproduct / scores1.length;
        result = cov_x_y / (pop_sd_x * pop_sd_y);
        return result;
    }

    private static void testPearson() throws IOException {

        String file = "/Users/jrobinso/IGV/hic/pearson/OE.txt";
        String ofile = "/Users/jrobinso/IGV/hic/pearson/out3.txt";
        parseOEFile(file, ofile);

//        To compute the Pearson's product-moment correlation between two double arrays x and y, use:
//
//        new PearsonsCorrelation().correlation(x, y)
//


    }

    private static void parseOEFile(String oeFile, String ofile) throws IOException {

        BufferedReader br = null;
        double[][] vectors = null;
        List<Integer> zeroRows = new ArrayList();
        try {
            br = new BufferedReader(new FileReader(oeFile));
            String nexLine;

            String header = br.readLine();
            String[] tokens = header.split("\t");
            int nCols = tokens.length - 1;


            int rowNumber = 0;
            vectors = new double[nCols][];
            while ((nexLine = br.readLine()) != null) {

                boolean allZeroes = true;
                double[] vector = new double[nCols];
                tokens = nexLine.split("\t");
                if (tokens.length < nCols + 1) {
                    continue;
                }
                for (int i = 1; i < nCols; i++) {
                    try {
                        vector[i - 1] = Double.parseDouble(tokens[i]);
                        if (vector[i - 1] != 0) {
                            allZeroes = false;
                        }
                    } catch (NumberFormatException e) {
                        vector[i - 1] = 0;
                    }
                }
                if (allZeroes) {
                    zeroRows.add(rowNumber);
                }
                vectors[rowNumber] = vector;
                rowNumber++;
            }
        } finally {
            if (br != null) br.close();
        }

        // remove allzero columns. This is a AxA matrix, rows == columns
        int nZeroes = zeroRows.size();
        int vecIdx = 0;
        for (double[] vec : vectors) {
            double[] newVec = new double[vec.length - nZeroes];
            int newIdx = 0;
            for (int i = 0; i < vec.length; i++) {
                if (!zeroRows.contains(i)) {
                    newVec[newIdx] = vec[i];
                    newIdx++;
                }
            }
            vectors[vecIdx] = newVec;
            vecIdx++;
        }

        computePearson(vectors, zeroRows, ofile);
    }

    private static void computePearson(double[][] inputVectors, List<Integer> zeroIndeces, String ofile) throws IOException {

        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(ofile)));

            for (int i = 0; i < inputVectors.length; i++) {
                if (zeroIndeces.contains(i)) {
                    for (int j = 0; j < inputVectors.length; j++) {
                        pw.print("0\t");
                    }
                } else {
                    double[] v1 = inputVectors[i];
                    for (int j = 0; j < inputVectors.length; j++) {
                        if (j == i) {
                            pw.print(1);
                        } else if (zeroIndeces.contains(j)) {
                            pw.print(0);
                        } else {
                            double[] v2 = inputVectors[j];
                            double p = getPearsonCorrelation(v1, v2);
                            pw.print(p);
                        }

                        pw.print("\t");
                    }
                }
                pw.println();
            }


        } finally {
            if (pw != null) pw.close();
        }

    }
}
