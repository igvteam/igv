package org.broad.igv.hic.tools;

import org.apache.batik.util.DoublyLinkedList;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.broad.igv.hic.MainWindow;
import org.broad.igv.hic.data.Chromosome;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.BufferedInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

/**
 * Computes an "expected" density vector.  Essentially there are 3 steps to using this class
 *
 * (1) instantiate it with a collection of Chromosomes (representing a genome) and a grid size
 * (2) loop through the pair data,  calling addDistance for each pair, to accumlate all counts
 * (3) when data loop is complete, call computeDensity to do the calculation
 *
 * Methods are provided to save the result of the calculation to a binary file, and restore it.  See the
 * DensityUtil class for example usage.
 *
 * @author Jim Robinson
 * @date 11/27/11
 */
public class DensityCalculation {


    private long totalLen ;  // Total length of the genome (sum of all chromosome lengths)
    private int gridSize;
    private int numberOfBins;
    private double[] actualDistances;
    private double[] possibleDistances;
    private double[] density;
    private double[] densityAvg;
    private Chromosome[] chromosomes;

    int totalCounts;

    // Map of chromosome index -> total count for that chromosome
    Map<Integer, Integer> chromosomeCounts;

    // Map of chromosome index -> "normalization factor", essentially a fudge factor to make
    // the "expected total"  == observed total
    private LinkedHashMap<Integer, Double> normalizationFactors;


    /**
     * Instantiate a DensityCalculation.  This constructor is used to compute the "expected" density from pair data.
     *
     * @param chromosomes
     * @param gridSize
     */
    DensityCalculation(Chromosome[] chromosomes, int gridSize) {
        this.gridSize = gridSize;
        this.chromosomes = chromosomes;

        totalLen = 0;
        for(Chromosome chromosome : chromosomes) {
            totalLen += chromosome.getSize();
        }

        numberOfBins = (int) (totalLen / gridSize) + 1;
        actualDistances = new double[numberOfBins];
        Arrays.fill(actualDistances, 0);
        chromosomeCounts = new HashMap();
        normalizationFactors = new LinkedHashMap<Integer, Double>();

    }


    /**
     * Read a previously calculaed density calculation from an input stream.
     *
     * @param is
     */
    public DensityCalculation(InputStream is) {
        LittleEndianInputStream les = new LittleEndianInputStream(new BufferedInputStream(is));
        try {
            read(les);
        } catch (IOException e) {
            System.err.println("Error reading density file");
            e.printStackTrace();
        }
    }


    /**
     * Add an observed distance.  This is called for each pair in the dataset
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

        // Compute "possible distances" for each bin.  I'm not sure I'm buying this but that's what the Python does.
        possibleDistances = new double[numberOfBins];
        for (Chromosome chr : chromosomes) {
            if (chr == null) continue;
            int nChrBins = chr.getSize() / gridSize;
            for (int i = 0; i < nChrBins; i++) {
                possibleDistances[i] += (nChrBins - i);
            }
        }

        // Compute the non-smoothed "density",  which is the actual count / possible count for each bin
        density = new double[numberOfBins];
        densityAvg = new double[numberOfBins];
        for (int i = 0; i < numberOfBins; i++) {
            density[i] = actualDistances[i] / possibleDistances[i];
            densityAvg[i] = density[i];  // <= initial value only,  this is "smoothed" below
        }


        // Smooth in 3 stages,  the window sizes are tuned to human.

        // Smooth (1)
        final int smoothingWidow1 = 15000000;
        int start = smoothingWidow1 / gridSize;
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



        // Compute fudge factors for each chromosome so the total "expected" count for that chromosome == the observed
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
                    // Skip.  Not sure why we should have these actually.

                } else {
                    expectedCount += (nGrids - n) * v;
                }
            }

            double observedCount = (double) chromosomeCounts.get(chr.getIndex());

            double f = expectedCount / observedCount;
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


    /**
     * Output the density calculation to a binary file.
     *
     * @param os
     * @throws IOException
     */
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

    /**
     * Read the contents of a previously saved calculation.
     *
     * @param is
     * @throws IOException
     */
    public void read(LittleEndianInputStream is) throws IOException {
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
        for (int i = 0; i < nDensities; i++) {
            densityAvg[i] = is.readDouble();
        }


    }


}


/*

--- Code above based on the following Python

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
