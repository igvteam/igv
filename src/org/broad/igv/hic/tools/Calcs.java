package org.broad.igv.hic.tools;

import org.broad.igv.hic.data.Chromosome;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author Jim Robinson
 * @date 11/27/11
 */
public class Calcs {

    int totalLen = 300000000;
    int gridSize = 1000000;  // 1MB
    int numberOfBins;
    double[] actualDistances;
    double[] possibleDistances;
    double[] density;
    double[] densityAvg;
    Chromosome[] chromosomes = HiCTools.b37Chromosomes;


    public static void main(String[] args) throws IOException {
        String path = "/Users/jrobinso/IGV/hic/formattedalignment.txt.gz";
        AsciiPairIterator iter = new AsciiPairIterator(path);

        Calcs calcs = new Calcs();
        while (iter.hasNext()) {
            AlignmentPair pair = iter.next();
            if (pair.getChr1().equals(pair.getChr2())) {
                int dist = Math.abs(pair.getPos1() - pair.getPos2());
                calcs.addDistance(dist);
            }
        }
        calcs.computeDensity();
        calcs.printDensities();
    }

    Calcs() {
        numberOfBins = (int) (totalLen / gridSize) + 1;
        actualDistances = new double[numberOfBins];
        Arrays.fill(actualDistances, 0);

    }

    public void addDistance(int dist) {
        int bin = dist / gridSize;
        actualDistances[bin]++;
    }

    public void computeDensity() {
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


    }


    private void printDensities() {

        int center = gridSize / 2;
        for (int i = 0; i < numberOfBins; i++) {
            System.out.println(center
                    + "\t" + actualDistances[i]
                    + "\t" + possibleDistances[i]
                    // + "\t" + density[i]
                    + "\t" + densityAvg[i]);
            center += gridSize;
        }

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
