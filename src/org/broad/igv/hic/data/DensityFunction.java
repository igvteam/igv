package org.broad.igv.hic.data;

import org.broad.igv.hic.tools.DensityCalculation;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 12/5/11
 */
public class DensityFunction {

    double[] density;
    int gridSize;
    private int nPoints;
    private double sum;
    Map<Integer, Double> normFactors;

    public DensityFunction(int gridSize, double[] densities, Map<Integer, Double> normFactors) {
        this.gridSize = gridSize;
        this.density = densities;
        this.nPoints = densities.length;
        this.normFactors = normFactors;
        this.sum = 0;
        for (int i=0; i<nPoints; i++)
            this.sum += density[i];
    }

    public DensityFunction(DensityCalculation calculation) {
        this(calculation.getGridSize(), calculation.getDensityAvg(), calculation.getNormalizationFactors());
    }

    public double getSum() {
        return this.sum;
    }

    public double getDensity(int chrIdx, int distance) {

       // double normFactor = normFactors.containsKey(chrIdx) ? normFactors.get(chrIdx) : 1.0;
      //  normFactor *= .8;
        double normFactor = 1;
        // Norm factor essentially present for backwards compatibility
        if (distance >= nPoints) {

            return density[nPoints - 1] / normFactor;
        } else {
            return density[distance] / normFactor;
        }
    }

}
