package org.broad.igv.hic.tools;

import org.broad.igv.feature.Chromosome;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.IOException;
import java.util.*;

/**
 * Computes an "expected" density vector.  Essentially there are 3 steps to using this class
 *
 * (1) instantiate it with a collection of Chromosomes (representing a genome) and a grid size
 * (2) loop through the pair data,  calling addDistance for each pair, to accumulate all counts
 * (3) when data loop is complete, call computeDensity to do the calculation
 *
 * Methods are provided to save the result of the calculation to a binary file, and restore it.  See the
 * DensityUtil class for example usage.
 *
 * @author Jim Robinson
 * @since 11/27/11
 */
public class DensityCalculation {

    private int gridSize;
    private int numberOfBins;
    private double[] actualDistances;    // genome wide
    private double[] possibleDistances;  // genome wide
    private double[] densityAvg;
    private List<Chromosome>  chromosomes;

    /** Map of chromosome index -> total count for that chromosome    */
    private Map<Integer, Integer> chromosomeCounts;

    /** Map of chromosome index -> "normalization factor", essentially a fudge factor to make
      * the "expected total"  == observed total */
    private LinkedHashMap<Integer, Double> normalizationFactors;

    /**
     * Instantiate a DensityCalculation.  This constructor is used to compute the "expected" density from pair data.
     *
     * @param chromosomes List of chromosomes, mainly used for size
     * @param gridSize    Grid size, used for binning appropriately
     */
    DensityCalculation(List<Chromosome> chromosomes, int gridSize) {

        this.gridSize = gridSize;

        this.chromosomes = chromosomes;
        long totalLen = 0;
        for(Chromosome chromosome : chromosomes) {
            if (chromosome != null)
                totalLen += chromosome.getLength();
        }

        numberOfBins = (int) (totalLen / gridSize) + 1;
        actualDistances = new double[numberOfBins];
        Arrays.fill(actualDistances, 0);
        chromosomeCounts = new HashMap<Integer,Integer>();
        normalizationFactors = new LinkedHashMap<Integer, Double>();
    }


    /**
     * Read a previously calculated density calculation from an input stream.
     *
     * @param les Stream to read from
     */
    public DensityCalculation(LittleEndianInputStream les) {
        try {
            read(les);
        } catch (IOException e) {
            System.err.println("Error reading density file");
            e.printStackTrace();
        }
    }


    /**
     * Add an observed distance.  This is called for each pair in the data set
     *
     * @param chr  Chromosome where observed, so can increment count
     * @param dist Distance observed
     */
    public void addDistance(Integer chr, int dist) {

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
     * Compute the "density" -- port of python function getDensityControls().
     * The density is a measure of the average distribution of counts genome-wide for a ligated molecule.
     * The density will decrease as distance from the center diagonal increases.
     * First compute "possible distances" for each bin.
     * "possible distances" provides a way to normalize the counts. Basically it's the number of
     * slots available in the diagonal.  The sum along the diagonal will then be the count at that distance,
     * an "expected" or average uniform density.
     */
    public void computeDensity() {
        double[] density;


        possibleDistances = new double[numberOfBins];

        for (Chromosome chr : chromosomes) {
            if (chr == null) continue;
            int nChrBins = chr.getLength() / gridSize;

            // this is not actually the true "possible", which would count everything, but is off by a factor
            // of binSize, essentially.  it makes the numbers more reasonable and is what the Python code does.
            for (int i = 0; i < nChrBins; i++) {
                possibleDistances[i] += (nChrBins - i);
            }
            // Lots of zeros at the end of this, because no contacts


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
        int window = (int) (5 * (2000000f / gridSize));
        if (window == 0) window = 1;
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
        window = (int)(20 * (2000000f / gridSize));
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


            int len = chr.getLength();
            int nGrids = len / gridSize + 1;
            double expectedCount = 0;
            for (int n = 0; n < nGrids; n++) {
                final double v = densityAvg[n];
                if (Double.isNaN(v)) {
                    System.err.println("Density was NaN, this shouldn't happen; possibly due to smoothing non-human genome");
                }
                else {
                    // this is the sum of the diagonal for this particular chromosome.
                    // the value in each bin is multiplied by the length of the diagonal to get expected count
                    expectedCount += (nGrids - n) * v;
                }
            }

            double observedCount = (double) chromosomeCounts.get(chr.getIndex());
            double f = expectedCount / observedCount;

            normalizationFactors.put(chr.getIndex(), f);
        }
    }

    /**
     * Debugging method that prints out everything.
     * @param chrIndex Chromosome to print.
     */
    private void printDensities(int chrIndex) {
        double norm = normalizationFactors.get(chrIndex);

        long center = gridSize / 2;
        System.out.println(gridSize);
        for (int i = 0; i < numberOfBins; i++) {

            System.out.println(center/1000000f + "M " +
                    actualDistances[i]
                    + "\t" + possibleDistances[i]
                    // + "\t" + density[i]
                    + "\t" + densityAvg[i]     + "\t" + norm
                    + "\t" + densityAvg[i] / norm);
            center += gridSize;
        }

    }

    /**
     * Accessor for grid size
     * @return The grid size
     */
    public int getGridSize() {
        return gridSize;
    }

    /**
     * Accessor for the normalization factors
     * @return  The normalization factors
     */
    public LinkedHashMap<Integer, Double> getNormalizationFactors() {
        return normalizationFactors;
    }

    /**
     * Accessor for the densities
     * @return The densities
     */
    public double[] getDensityAvg() {
        return densityAvg;
    }

    /**
     * Output the density calculation to a binary file.
     *
     * @param os      Stream to output to
     * @throws IOException    If error while writing
     */
    public void outputBinary(LittleEndianOutputStream os) throws IOException {

        os.writeInt(gridSize);

        int nonNullChrCount = 0;
        for (Chromosome chr : chromosomes) {
            if (chr != null) nonNullChrCount++;
        }
        os.writeInt(nonNullChrCount);

        // Chromosome indices
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
            os.writeDouble(normFactor == null ? 0 : normFactor);
        }

        // Expected value (densityAvg)
        os.writeInt(densityAvg.length);
        for (double aDensityAvg : densityAvg) {
            os.writeDouble(aDensityAvg);
        }
    }

    /**
     * Read the contents of a previously saved calculation.
     *
     * @param is  Stream to read
     * @throws IOException      If error while reading stream
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
        normalizationFactors = new LinkedHashMap<Integer,Double>(nChromosomes);
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
