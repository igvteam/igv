package org.broad.igv.hic.data;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.tools.ExpectedValueCalculation;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Utility holder for Density calculation, for O/E maps.
 *
 * @author Jim Robinson
 * @author Neva Cherniavsky
 * @since 8/27/12
 */
public class DensityFunction {

    int binSize;

    String unit;

    public Map<Integer, Double> normFactors;

    double[] density;

    public DensityFunction(String unit, int binSize, double[] density, Map<Integer, Double> normFactors) {
        this.unit = unit;
        this.binSize = binSize;
        this.normFactors = normFactors;
        this.density = density;
    }

    /**
     * Gets the expected value, distance and coverage normalized
     *
     * @param chrIdx   Chromosome index
     * @param distance Distance from diagonal in bins
     * @return Expected value, distance and coverage normalized
     */
    public double getDensity(int chrIdx, int distance) {

        double normFactor = 1.0;
        if (normFactors != null && normFactors.containsKey(chrIdx)) {
            normFactor = normFactors.get(chrIdx);
        }

        if (distance >= density.length) {

            return density[density.length - 1] / normFactor;
        } else {
            return density[distance] / normFactor;
        }
    }


}
