package org.broad.igv.hic.tools;

import org.apache.batik.util.DoublyLinkedList;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.broad.igv.hic.MainWindow;
import org.broad.igv.hic.data.Chromosome;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.DataOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * @author Jim Robinson
 * @date 11/27/11
 */
public class DensityCalculation {


    private int totalLen = 300000000;
    private int gridSize = 1000000;  // 1MB
    private int numberOfBins;
    private double[] actualDistances;
    private double[] possibleDistances;
    private double[] density;
    private double[] densityAvg;
    private Chromosome[] chromosomes;

    int totalCounts;

    // chr -> count
    Map<Integer, Integer> chromosomeCounts;

    // chr -> norm
    private LinkedHashMap<Integer, Double> normalizationFactors;

    public DensityCalculation() {
    }

    DensityCalculation(Chromosome[] chromosomes, int gridSize) {
        this.gridSize = gridSize;
        this.chromosomes = chromosomes;
        numberOfBins = (totalLen / gridSize) + 1;
        actualDistances = new double[numberOfBins];
        Arrays.fill(actualDistances, 0);
        chromosomeCounts = new HashMap();
        normalizationFactors = new LinkedHashMap<Integer, Double>();

    }


    public void outputBinary(LittleEndianOutputStream os) throws IOException {

        os.writeInt(gridSize);

        int nonNullChrCount = 0;
        for (Chromosome chr : chromosomes) {
            if (chr != null) nonNullChrCount++;
        }
        os.writeInt(nonNullChrCount);

        // Chromosome indeces
        for (Chromosome chr : chromosomes) {
            if (chr == null) continue;
            os.writeInt(chr.getIndex());
        }

        // Normalization factors
        for (Chromosome chr : chromosomes) {
            if (chr == null) continue;
            Integer idx = chr.getIndex();
            Double normFactor = normalizationFactors.get(idx);
            os.writeInt(idx);
            os.writeDouble(normFactor == null ? 0 : normFactor.doubleValue());
        }

        // Expected value (densityAvg)
        os.writeInt(densityAvg.length);
        for (int i = 0; i < densityAvg.length; i++) {
            os.writeDouble(densityAvg[i]);
        }
    }

    public void fill(LittleEndianInputStream is) throws IOException {
        gridSize = is.readInt();
        int nChromosomes = is.readInt();

        // Chromosome indexes
        Integer[] chrIndexes = new Integer[nChromosomes];
        for (int i = 0; i < nChromosomes; i++) {
            chrIndexes[i] = is.readInt();
        }

        // Normalization factors
        normalizationFactors = new LinkedHashMap(nChromosomes);
        for (int i = 0; i < nChromosomes; i++) {
            Integer chrIdx = is.readInt();
            double normFactor = is.readDouble();
            normalizationFactors.put(chrIdx, normFactor);
        }

        // Densities
        int nDensities = is.readInt();
        densityAvg = new double[nDensities];
        for(int i=0; i<nDensities; i++) {
            densityAvg[i] = is.readDouble();
        }


    }

    /**
     * Add an observed distance
     *
     * @param chr
     * @param dist
     */
    public void addDistance(Integer chr, int dist) {

        totalCounts++;

        Integer count = chromosomeCounts.get(chr);
        if (count == null) {
            chromosomeCounts.put(chr, 1);
        } else {
            chromosomeCounts.put(chr, count + 1);
        }
        int bin = dist / gridSize;
        actualDistances[bin]++;
    }

    /**
     * Compute the "density" -- port of python function getDensityControls()
     */
    public void computeDensity() {

        // Compute possible distances
        possibleDistances = new double[numberOfBins];
        for (Chromosome chr : chromosomes) {
            if (chr == null) continue;
            int nChrBins = chr.getSize() / gridSize;
            for (int i = 0; i < nChrBins; i++) {
                possibleDistances[i] += (nChrBins - i);
            }
        }

        density = new double[numberOfBins];
        densityAvg = new double[numberOfBins];
        for (int i = 0; i < numberOfBins; i++) {
            density[i] = actualDistances[i] / possibleDistances[i];
            densityAvg[i] = density[i];
        }

        // Smooth (1)
        int start = 15000000 / gridSize;
        int window = 5 * (2000000 / gridSize);
        for (int i = start; i < numberOfBins; i++) {
            int kMin = i - window;
            int kMax = Math.min(i + window, numberOfBins);
            double sum = 0;
            for (int k = kMin; k < kMax; k++) {
                sum += density[k];
            }
            densityAvg[i] = sum / (kMax - kMin);
        }

        // Smooth (2)
        start = 70000000 / gridSize;
        window = 20 * (2000000 / gridSize);
        for (int i = start; i < numberOfBins; i++) {
            int kMin = i - window;
            int kMax = Math.min(i + window, numberOfBins);
            double sum = 0;
            for (int k = kMin; k < kMax; k++) {
                sum += density[k];
            }
            densityAvg[i] = sum / (kMax - kMin);
        }

        // Smooth (3)
        start = 170000000 / gridSize;
        for (int i = start; i < numberOfBins; i++) {
            densityAvg[i] = densityAvg[start];
        }


        //4 3 2 1
        //  4 3 2
        //    4 3
        //      4
        // g0 = 4,  g1 = 3,  g2 = 2,  g3 = 1,  N=4
        // total = N*g0 + (N-1)*g1 + (N-2)*g2 + (N-3)*g3

        for (Chromosome chr : chromosomes) {

            if (chr == null || !chromosomeCounts.containsKey(chr.getIndex())) {
                continue;
            }

            int len = chr.getSize();
            int nGrids = len / gridSize + 1;
            double expectedCount = 0;
            for (int n = 0; n < nGrids; n++) {
                final double v = densityAvg[n];
                if (Double.isNaN(v)) {
                    System.out.println(n);
                } else {
                    expectedCount += (nGrids - n) * v;
                }
            }

            double observedCount = (double) chromosomeCounts.get(chr.getIndex());

            double f = observedCount / expectedCount;
            System.out.println(chr.getName() + "\t" + f);

            normalizationFactors.put(chr.getIndex(), f);
        }
    }


    private void printDensities(int chrIndex) {


        double norm = normalizationFactors.get(chrIndex);

        int center = gridSize / 2;
        for (int i = 0; i < numberOfBins; i++) {
            System.out.println(center
                    + "\t" + actualDistances[i]
                    + "\t" + possibleDistances[i]
                    // + "\t" + density[i]
                    + "\t" + densityAvg[i]
                    + "\t" + densityAvg[i] * norm);
            center += gridSize;
        }

    }


    public int getGridSize() {
        return gridSize;
    }

    public LinkedHashMap<Integer, Double> getNormalizationFactors() {
        return normalizationFactors;
    }

    public int getNumberOfBins() {
        return numberOfBins;
    }

    public double[] getDensityAvg() {
        return densityAvg;
    }


}


/*
gridSize => grid (or bin) size  == 10^6
actualDistances => array of actual distances,  each element represents a bin
possibleDistances => array of possible distances, each element represents a bin
locs => chromosome lengths
jdists => outer distances between pairs

for each jdist
  actualDistance[jdist]++;


for each chromosome
  chrlen = chromosome length
  numberOfBins = chrlen / gridSize
  for each i from 0 to numberOfBins
     possibleDistances[i] += (numberOfBins - i)


for each i from 0 to maxGrid
  density[i] = actualDistance[i] / possibleDistances[i]


for each i from 0 to len(density)
 density_avg[i] = density[i]

for each i from 15000000/gridsize  to  len(density_avg)
  sum1 = 0
  for each k from (i - 5*((2*10^6) / gridSize)  to  (i + 5*((2*10^6)/gridsize))
     sum1 += density[k]
  density_avg[i] = sum1 / (10*((2*10^6)/gridsize))

for each i from 70000000/gridsize  to  len(density_avg)
  sum2 = 0
  for each k from (i - 20*((2*10^6) / gridSize)  to  (i + 20*((2*10^6)/gridsize))
     sum2 += density[k]
  density_avg[i] = sum2 / (40*((2*10^6)/gridsize))

for each i from 170000000/gridsize  to  len(density_avg)
  density_avg[i]=density_avg[170000000/gridsize]

*/
