/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tdf.Accumulator;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.util.collections.IntArrayList;

import java.util.*;

/**
 * Summarize (using a windowing function) numeric data points which are associated
 * with locations on a genome. Stored by chromosome
 * @author jrobinso
 */
public class GenomeSummaryData {

    private static Logger log = Logger.getLogger(GenomeSummaryData.class);

    // Genome coordinates are in kilobases
    private static final double locationUnit = 1000.0;

    /**
     * Number of virtual pixels
     */
    int nPixels = 1000;

    Genome genome;

    String[] samples;

    /**
     * Chromosome name -> [sample name -> list of data values]
     */
    Map<String, Map<String, FloatArrayList>> dataMap = new HashMap<String, Map<String, FloatArrayList>>();

    /**
     * Chromosome name -> list of start locations
     */
    Map<String, IntArrayList> locationMap;

    /**
     * start locations of currently relevant ordered list of chromosomes, spanning whole genomes
     */
    int[] locations;


    /**
     * sample name -> list of sample data
     */
    Map<String, float[]> data;

    int nDataPts = 0;

    Set<String> skippedChromosomes = new HashSet<String>();


    /**
     * Scale in kilobases / pixel
     */
    double scale;

    public GenomeSummaryData(Genome genome, String[] samples) {
        this.genome = genome;
        this.samples = samples;
        scale = (genome.getNominalLength() / locationUnit) / nPixels;

        List<String> chrNames = genome.getLongChromosomeNames();
        locationMap = new HashMap<String, IntArrayList>();
        dataMap = new HashMap<String, Map<String, FloatArrayList>>();
        for (String chr : chrNames) {
            locationMap.put(chr, new IntArrayList(nPixels / 10));
            dataMap.put(chr, new HashMap<String, FloatArrayList>());
            for (String s : samples) {
                dataMap.get(chr).put(s, new FloatArrayList(nPixels / 10));
            }
        }
    }

    /**
     * Changes scale of summary, ie zoom in or out
     * Mainly for testing, can't use after adding any data
     * @param scale
     */
    void setScale(double scale){
        if(nDataPts > 0) throw new IllegalStateException("Can't alter scale after adding data");
        this.scale = scale;
        nPixels = (int) (((double) this.genome.getNominalLength() / locationUnit) / scale);
    }


    /**
     * Add data to be condensed for the whole genome view
     *
     * @param chr
     * @param locs Genomic positions
     * @param sampleData  Map of sample name -> array of values.
     */
    public void addData(String chr, int[] locs, Map<String, float[]> sampleData) {

        IntArrayList locations = locationMap.get(chr);
        if (locations == null) {
            if (!skippedChromosomes.contains(chr)) {
                skippedChromosomes.add(chr);
                log.info("Skipping data for: " + chr);
            }
            return;
        }

        int lastPixel = -1;
        int lastGenomeLocation = -1;
        Map<String, Accumulator> dataPoints = new HashMap<String, Accumulator>();

        for (int i = 0; i < locs.length; i++) {

            int genomeLocation = genome.getGenomeCoordinate(chr, locs[i]);
            int pixel = (int) (genomeLocation / scale);
            if (lastPixel >= 0 && pixel != lastPixel) {
                locations.add(lastGenomeLocation);
                finishLastLocation(chr, dataPoints);
            }

            for (String s : samples) {
                float[] data = sampleData.get(s);
                Accumulator dp = dataPoints.get(s);
                if (dp == null) {
                    dp = new Accumulator(WindowFunction.mean);
                    dataPoints.put(s, dp);
                }
                try {
                    dp.add(1, data[i], null);
                } catch (Exception e) {
                    log.error("Error adding to GenomeSummaryData", e);
                }
            }

            lastPixel = pixel;
            lastGenomeLocation = genomeLocation;
        }

        locations.add(lastGenomeLocation);
        finishLastLocation(chr, dataPoints);
    }

    /**
     * Mark the previous genomic location as having been completely summarized
     * @param chr
     * @param dataPoints  Map sample -> accumulator, which stored data temporarily being accumulated at a given genome location
     */
    private void finishLastLocation(String chr, Map<String, Accumulator> dataPoints) {
        nDataPts++;

        for (String s : dataMap.get(chr).keySet()) {
            Accumulator dp = dataPoints.get(s);
            dp.finish();
            dataMap.get(chr).get(s).add(dp.getValue());
        }
        dataPoints.clear();
    }

    public int[] getLocations() {
        if (locations == null) {
            createDataArrays();
        }

        return locations;
    }


    public float[] getData(String sample) {
        if (!data.containsKey(sample)) {
            createDataArrays();
        }
        return data.get(sample);

    }

    /**
     * Recalculate:
     * 0. Start locations for plotting. Shared across all samples
     * 1. Summary data for a given sample, across all stored chromosomes
     */
    private synchronized void createDataArrays() {
        locations = new int[nDataPts];
        int offset = 0;
        List<String> chrNames = genome.getLongChromosomeNames();
        for (String chr : chrNames) {
            int[] chrLocs = locationMap.get(chr).toArray();
            System.arraycopy(chrLocs, 0, locations, offset, chrLocs.length);
            offset += chrLocs.length;
        }

        data = new HashMap<String, float[]>();
        for (String s : samples) {
            float[] sampleData = new float[nDataPts];
            offset = 0;
            for (String chr : chrNames) {
                float[] chrData = dataMap.get(chr).get(s).toArray();
                System.arraycopy(chrData, 0, sampleData, offset, chrData.length);
                offset += chrData.length;
            }
            data.put(s, sampleData);
        }

        locationMap.clear();
        dataMap.clear();
    }

}
