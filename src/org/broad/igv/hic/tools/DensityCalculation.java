package org.broad.igv.hic.tools;

import org.broad.igv.feature.Chromosome;
import org.broad.tribble.util.LittleEndianInputStream;
import java.io.IOException;
import java.util.*;

/**
 * Computes an "expected" density vector.  Essentially there are 3 steps to using this class
 *
 * (1) instantiate it with a collection of Chromosomes (representing a genome) and a grid size
 * (2) loop through the pair data,  calling addDistance for each pair, to accumulate all counts
 * (3) when data loop is complete, call computeDensity to do the calculation
 *     - with new version, this will also compute coverage normalization
 *
 * Methods are provided to save the result of the calculation to a binary file, and restore it.  See the
 * DensityUtil class for example usage.
 *
 * @author Jim Robinson
 * @since 11/27/11
 */
public class DensityCalculation {

    private int      gridSize;
    private int      numberOfBins;
    private long     totalReads;
    private boolean  isNewVersion;
    /** Genome wide count of binned reads at a given distance */
    private double[] actualDistances     = null;
    /** Genome wide binned possible distances */
    private double[] possibleDistances   = null;
    /** Sum of counts in a row, used for coverage normalization */
    private double[] rowSums             = null;
    /** Coverage normalizations for each row (= column) */
    private double[] coverageNorms       = null;
    /** Expected count at a given binned distance from diagonal */
    private double[] densityAvg          = null;
    /** Chromosome in this genome, needed for normalizations */
    private List<Chromosome> chromosomes = null;

    /** Map of chromosome index -> total count for that chromosome    */
    private Map<Integer, Integer> chromosomeCounts;

    /** Map of chromosome index -> "normalization factor", essentially a fudge factor to make
      * the "expected total"  == observed total */
    private LinkedHashMap<Integer, Double> normalizationFactors;

    /** Stores restriction site fragment information for fragment maps */
    private FragmentCalculation fragmentCalculation;

    /**
     * Instantiate a DensityCalculation.  This constructor is used to compute the "expected" density from pair data.
     *
     * @param chromosomes List of chromosomes, mainly used for size
     * @param gridSize    Grid size, used for binning appropriately
     */
    DensityCalculation(List<Chromosome> chromosomes, int gridSize, FragmentCalculation fragmentCalculation, boolean isNewVersion) {

        this.gridSize = gridSize;

        this.chromosomes = chromosomes;
        long totalLen = 0;
        for (Chromosome chromosome : chromosomes) {
            if (chromosome != null) {
                if (gridSize == 1)
                    totalLen += fragmentCalculation.getNumberFragments(chromosome);
                else
                    totalLen += chromosome.getLength();
            }
        }
        if (gridSize == 1)
            numberOfBins = (int) totalLen;
        else
            numberOfBins = (int) (totalLen / gridSize) + 1;
        actualDistances = new double[numberOfBins];
        rowSums = new double[numberOfBins];
        coverageNorms = new double[numberOfBins];
        Arrays.fill(coverageNorms, 1);
        Arrays.fill(actualDistances, 0);
        Arrays.fill(rowSums, 0);
        chromosomeCounts = new HashMap<Integer,Integer>();
        normalizationFactors = new LinkedHashMap<Integer, Double>();
        totalReads = 0;
        this.isNewVersion = isNewVersion;
        this.fragmentCalculation = fragmentCalculation;
    }

    /**
     * Read a previously calculated density calculation from an input stream.
     *
     * @param les Stream to read from
     * @param isNewVersion New version stores the coverage normalization as well
     */
    public DensityCalculation(LittleEndianInputStream les, boolean isNewVersion) {
        try {
            read(les, isNewVersion);
            this.isNewVersion = isNewVersion;
        } catch (IOException e) {
            System.err.println("Error reading density file");
            e.printStackTrace();
        }
    }

    /**
     * Set list of chromosomes; need to do this when reading from file
     * @param chromosomes1    Array of chromosomes to set
     */
    public void setChromosomes(Chromosome[] chromosomes1) {
        chromosomes = Arrays.asList(chromosomes1);
    }

    /**
     * Add an observed distance.  This is called for each pair in the data set
     *
     * @param chr  Chromosome where observed, so can increment count
     * @param pos1 Position1 observed
     * @param pos2 Position2 observed
     */
    public void addDistance(Integer chr, int pos1, int pos2) {
        int dist;
        Integer count = chromosomeCounts.get(chr);
        if (count == null) {
            chromosomeCounts.put(chr, 1);
        } else {
            chromosomeCounts.put(chr, count + 1);
        }
        if (gridSize == 1) {
            int bin1 = fragmentCalculation.getBin(chr, pos1);
            int bin2 = fragmentCalculation.getBin(chr, pos2);
            dist = Math.abs(bin1 - bin2);
        }
        else {
            dist = Math.abs(pos1 - pos2);
        }
        int bin = dist / gridSize;

        if (isNewVersion) {

            int bin1 = getGenomicRowBin(chr, pos1);
            int bin2 = getGenomicRowBin(chr, pos2);

            actualDistances[bin]+= 1/(coverageNorms[bin1]*coverageNorms[bin2]);
        }
        else
            actualDistances[bin]++;

    }

    /**
     * Returns normalized version of count, based on place in genome-wide matrix
     * @param count  Sum to normalize
     * @param chr1   First chromosome
     * @param pos1   First position
     * @param chr2   Second chromosome
     * @param pos2   Second position
     * @return  If old version and we don't have normalization, return count.  Else look at genome-wide bin
     * and divide by the normalizations for the row and column.
     */
    public double getNormalizedCount(int count, int chr1, int pos1, int chr2, int pos2) {
        if (coverageNorms == null) {
            return count;
        }
        int bin1 = getGenomicRowBin(chr1,pos1);
        int bin2 = getGenomicRowBin(chr2,pos2);
        return (double)count/(coverageNorms[bin1]*coverageNorms[bin2]);
    }

    /**
     * Add to rowSums array, so that we can do coverage normalization later
     * @param pair Pair to add
     */
    public void addToRow(AlignmentPair pair) {

        int bin1 = getGenomicRowBin(pair.getChr1(), pair.getPos1());
        int bin2 = getGenomicRowBin(pair.getChr2(), pair.getPos2());

        rowSums[bin1]++;
        // don't double count on the diagonal
        if (bin1 != bin2)
            rowSums[bin2]++;

        totalReads++;
    }

    /**
     * Find row bin in big genomic matrix
     * @param chr Chromosome
     * @param pos Position
     * @return  Row bin in genomic matrix
     */
    public int getGenomicRowBin(int chr, int pos) {
        long rowsBefore = 0;

        for (int i = 1; i < chr; i++) {
            rowsBefore += chromosomes.get(i).getLength();
        }

        rowsBefore += pos;
        return (int)(rowsBefore / gridSize);
    }

    public void computeCoverageNormalization() {
        double rowMean = totalReads / rowSums.length;
        for (int i=0; i < rowSums.length; i++) {
            coverageNorms[i] = rowSums[i] / rowMean;
        }
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
        int trueNumBins = 0;

        possibleDistances = new double[numberOfBins];

        for (Chromosome chr : chromosomes) {
            Integer count = chromosomeCounts.get(chr.getIndex());
            // didn't see anything at all from a chromosome, then don't include it in possDists.
            if (count == null) continue;
            if (chr == null) continue;
            int nChrBins;
            if (gridSize == 1)
                nChrBins = fragmentCalculation.getNumberFragments(chr);
            else
                nChrBins = chr.getLength();
            nChrBins = nChrBins / gridSize;

            for (int i = 0; i < nChrBins; i++) {
                possibleDistances[i] += (nChrBins - i);

                if (i > trueNumBins)
                    trueNumBins = i;
            }


        }

        density = new double[trueNumBins];

        densityAvg = new double[trueNumBins];
        for (int i = 0; i < trueNumBins; i++) {
            density[i] = actualDistances[i] / possibleDistances[i];
            if (actualDistances[i] < 400) {
                double tmp = actualDistances[i];
                double poss = 0;
                int window = 0;
                while (tmp < 400) {
                    window++;
                    tmp = 0;
                    for (int j=i-window; j <= i+window; j++)
                        tmp += actualDistances[j];
                }
                tmp = 0;
                for (int j=i-window; j <= i+window; j++) {
                    tmp += actualDistances[j];
                    poss += possibleDistances[j];
                }
                densityAvg[i] = tmp/poss;
            }
            else
                densityAvg[i] = density[i];

        }

        // Compute fudge factors for each chromosome so the total "expected" count for that chromosome == the observed
        for (Chromosome chr : chromosomes) {

            if (chr == null || !chromosomeCounts.containsKey(chr.getIndex())) {
                continue;
            }
            int nGrids;
            if (gridSize == 1)
                nGrids = fragmentCalculation.getNumberFragments(chr);
            else
                nGrids = chr.getLength() / gridSize + 1;


            double expectedCount = 0;
            // ARGH THIS IS WRONG
            for (int n = 0; n < trueNumBins; n++) {
                final double v = densityAvg[n];
                if (Double.isNaN(v)) {
                    System.err.println("Density was NaN, this shouldn't happen");
                }
                else {
                    // this is the sum of the diagonal for this particular chromosome.
                    // the value in each bin is multiplied by the length of the diagonal to get expected count
                    if (nGrids > n)
                        expectedCount += (nGrids - n) * v;
                }
            }

            double observedCount = (double) chromosomeCounts.get(chr.getIndex());
            double f = expectedCount / observedCount;
            normalizationFactors.put(chr.getIndex(), f);
        }
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
     * Output the density calculation to a binary file buffer.
     *
     * @param buffer      Stream to output to
     * @throws IOException    If error while writing
     */
    public void outputBinary(Preprocessor.BufferedByteWriter buffer, boolean isNewVersion) throws IOException {

        buffer.putInt(gridSize);

        int nonNullChrCount = 0;
        for (Chromosome chr : chromosomes) {
            if (chr != null) nonNullChrCount++;
        }
        buffer.putInt(nonNullChrCount);

        // Chromosome indices
        for (Chromosome chr : chromosomes) {
            if (chr == null) continue;
            buffer.putInt(chr.getIndex());
        }

        // Normalization factors
        for (Chromosome chr : chromosomes) {
            if (chr == null) continue;
            Integer idx = chr.getIndex();
            Double normFactor = normalizationFactors.get(idx);
            buffer.putInt(idx);
            buffer.putDouble(normFactor == null ? 0 : normFactor);
        }

        // Expected value (densityAvg)
        buffer.putInt(densityAvg.length);
        for (double aDensityAvg : densityAvg) {
            buffer.putDouble(aDensityAvg);
        }
        // Coverage normalizations
        if (isNewVersion) {
            buffer.putInt(coverageNorms.length);
            for (double aRowNorm : coverageNorms) {
                buffer.putDouble(aRowNorm);
            }
        }
    }

    /**
     * Read the contents of a previously saved calculation.
     *
     * @param is  Stream to read
     * @throws IOException      If error while reading stream
     */
    public void read(LittleEndianInputStream is, boolean isNewVersion) throws IOException {
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

        if (isNewVersion) {
            // Norms
            int nNorms = is.readInt();
            coverageNorms = new double[nNorms];
            for (int i=0; i < nNorms; i++) {
                coverageNorms[i] = is.readDouble();

            }
        }

    }

}






// Smooth in 3 stages,  the window sizes are tuned to human.

//        // Smooth (1)
//        final int smoothingWidow1 = 15000000;
//        int start = smoothingWidow1 / gridSize;
//        int window = (int) (5 * (2000000f / gridSize));
//        if (window == 0) window = 1;
//        for (int i = start; i < numberOfBins; i++) {
//            int kMin = i - window;
//            int kMax = Math.min(i + window, numberOfBins);
//            double sum = 0;
//            for (int k = kMin; k < kMax; k++) {
//                sum += density[k];
//            }
//            densityAvg[i] = sum / (kMax - kMin);
//        }
//
//        // Smooth (2)
//        start = 70000000 / gridSize;
//        window = (int)(20 * (2000000f / gridSize));
//        for (int i = start; i < numberOfBins; i++) {
//            int kMin = i - window;
//            int kMax = Math.min(i + window, numberOfBins);
//            double sum = 0;
//            for (int k = kMin; k < kMax; k++) {
//                sum += density[k];
//            }
//            densityAvg[i] = sum / (kMax - kMin);
//        }
//
//        // Smooth (3)
//        start = 170000000 / gridSize;
//        for (int i = start; i < numberOfBins; i++) {
//            densityAvg[i] = densityAvg[start];
//        }


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
