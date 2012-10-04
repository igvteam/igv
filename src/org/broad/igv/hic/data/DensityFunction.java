package org.broad.igv.hic.data;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.tools.ExpectedValueCalculation;

import java.util.Map;

/**
 * Utility holder for Density calculation, for O/E maps.
 *
 * @author Jim Robinson
 * @author Neva Cherniavsky
 * @since 8/27/12
 */
public class DensityFunction {

    private ExpectedValueCalculation densityCalculation;

    /**
     * Constructor sets the density calculation.
     *
     * @param calculation Density calculation, containing distance expectation and coverage normalization
     */
    public DensityFunction(ExpectedValueCalculation calculation) {
        this.densityCalculation = calculation;
    }

    public void setChromosomes(Chromosome[] chromosomes) {
        densityCalculation.setChromosomes(chromosomes);
    }

    /**
     * Gets the expected value, distance and coverage normalized
     * @param chrIdx Chromosome index
     * @param distance Distance from diagonal
     * @return  Expected value, distance and coverage normalized
     */
    public double getDensity(int chrIdx, int distance) {
       Map<Integer, Double> normFactors = densityCalculation.getNormalizationFactors();
       double normFactor = normFactors.containsKey(chrIdx) ? normFactors.get(chrIdx) : 1.0;
       normFactor = 1; // change this in the future but right now these are really messed up.
       double density[] = densityCalculation.getDensityAvg();
        if (distance >= density.length) {

            return density[density.length - 1] / normFactor;
        } else {
            return density[distance] / normFactor;
        }
    }


}
