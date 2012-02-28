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


    //String densityFile = "HumanAugust.densities.txt";
    int[] positions;
    double[] density;
    int gridSize;
    private int nPoints;
    Map<Integer, Double> normFactors;

    public DensityFunction(int gridSize, double[] densities, Map<Integer, Double> normFactors) {
        this.gridSize = gridSize;
        this.density = densities;
        this.nPoints = densities.length;
        this.normFactors = normFactors;
    }

    public DensityFunction(DensityCalculation calculation) {
        this(calculation.getGridSize(), calculation.getDensityAvg(), calculation.getNormalizationFactors());
    }
    
    public double getDensity(int chrIdx, int distance) {

        double normFactor = normFactors.containsKey(chrIdx) ? normFactors.get(chrIdx) : 1.0;
        int grid = distance;
        if (grid >= nPoints) {

            return density[nPoints - 1] / normFactor;
        } else {
            return density[grid] / normFactor;
        }
    }

}
