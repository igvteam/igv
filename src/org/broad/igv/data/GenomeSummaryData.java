/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
 * @author jrobinso
 */
public class GenomeSummaryData {

    private static Logger log = Logger.getLogger(GenomeSummaryData.class);

    // Genome coordinates are in KB
    private static final double locationUnit = 1000.0;

    /**
     * Number of virtual pixels
     */
    int nPixels = 1000;

    Genome genome;

    String[] samples;

    Map<String, Map<String, FloatArrayList>> dataMap = new HashMap();

    Map<String, IntArrayList> locationMap;

    int[] locations;

    Map<String, float[]> data;

    int nDataPts = 0;

    Set<String> skippedChromosomes = new HashSet();


    /**
     * Scale in KB / pixel
     */
    double scale;

    public GenomeSummaryData(Genome genome, String[] samples) {
        this.genome = genome;
        this.samples = samples;
        scale = (genome.getNominalLength() / locationUnit) / nPixels;

        List<String> chrNames = genome.getLongChromosomeNames();
        locationMap = new HashMap();
        dataMap = new HashMap();
        for (String chr : chrNames) {
            locationMap.put(chr, new IntArrayList(nPixels / 10));
            dataMap.put(chr, new HashMap());
            for (String s : samples) {
                dataMap.get(chr).put(s, new FloatArrayList(nPixels / 10));
            }
        }
    }


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
        Map<String, Accumulator> dataPoints = new HashMap();

        for (int i = 0; i < locs.length; i++) {

            int genomeLocation = genome.getGenomeCoordinate(chr, locs[i]);
            int pixel = (int) (genomeLocation / scale);
            if (i > 0 && pixel != lastPixel) {
                nDataPts++;

                locations.add(genomeLocation);
                for (String s : dataMap.get(chr).keySet()) {
                    Accumulator dp = dataPoints.get(s);
                    dp.finish();
                    dataMap.get(chr).get(s).add(dp.getValue());
                }
                dataPoints.clear();
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
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

            lastPixel = pixel;
        }
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

    private synchronized void createDataArrays() {
        locations = new int[nDataPts];
        int offset = 0;
        List<String> chrNames = genome.getLongChromosomeNames();
        for (String chr : chrNames) {
            int[] chrLocs = locationMap.get(chr).toArray();
            System.arraycopy(chrLocs, 0, locations, offset, chrLocs.length);
            offset += chrLocs.length;
        }

        data = new HashMap();
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
