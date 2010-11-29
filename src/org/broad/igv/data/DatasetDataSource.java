/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;

import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class DatasetDataSource extends AbstractDataSource {

    private Logger log = Logger.getLogger(DatasetDataSource.class);

    String trackId;
    Dataset dataset;
    GenomeSummaryData genomeSummaryData;

    /**
     * Constructs ...
     *
     * @param genome
     * @param trackId
     * @param dataset
     */
    public DatasetDataSource(Genome genome, String trackId, Dataset dataset) {
        this.trackId = trackId;
        this.dataset = dataset;

        // TODO -- remove this "instanceof" hack
        if (genome.getHomeChromosome().equals(Globals.CHR_ALL)) {
            if (dataset instanceof IGVDataset) {
                genomeSummaryData = ((IGVDataset) dataset).getGenomeSummary();
            } else {
                genomeSummaryData = new GenomeSummaryData(genome, new String[]{trackId});
                for (Chromosome chr : genome.getChromosomes()) {
                    int[] startLocations = dataset.getStartLocations(chr.getName());
                    if (!chr.getName().equals(Globals.CHR_ALL) && (startLocations != null) && (startLocations.length > 0)) {
                        Map<String, float[]> dMap = new HashMap();
                        dMap.put(trackId, dataset.getData(trackId, chr.getName()));
                        genomeSummaryData.addData(chr.getName(), startLocations, dMap);
                    }
                }
            }
        }

        //if (genomeSummaryData.getLocations() == null || genomeSummaryData.getLocations().length == 0) {
        //    throw new RuntimeException("None of the chromosomes in the dataset were found on the loaded genome"
        //            + "(" + genome.getName() + ").");
        //}
    }

    @Override
    protected int getNumZoomLevels(String chr) {
        return 0;
    }

    @Override
    protected DataTile getRawData(String chr, int startLocation, int endLocation) {
        int[] startLocs = null;
        int[] endLocs = null;
        float[] data = null;
        String[] features = null;

        if (chr.equals(Globals.CHR_ALL) && genomeSummaryData != null) {
            startLocs = genomeSummaryData.getLocations();
            endLocs = startLocs;
            data = genomeSummaryData.getData(trackId);
        } else {
            startLocs = dataset.getStartLocations(chr);
            endLocs = dataset.getEndLocations(chr);
            data = dataset.getData(trackId, chr);
            features = dataset.getFeatureNames(chr);
        }
        if (startLocs == null) {
            return null;
        }

        assert (startLocs.length == endLocs.length);
        assert (startLocs.length == data.length);

        int maxIndex = startLocs.length - 1;

        // Until we can guarantee sorted values just use the index limits
        int startIdx = 0;    // Math.min(maxIndex,  DataUtils.getIndexBefore(startLocs, startLocation));
        int endIdx = maxIndex;    // Math.min(maxIndex, DataUtils.getIndexBefore(startLocs, endLocation) + 1);
        int[] tmpStart = startLocs;
        int[] tmpEnd = endLocs;
        float[] tmpData = data;
        String[] tmpFeatures = features;
        int nPts = endIdx - startIdx + 1;

        if ((tmpStart == null) || (tmpData == null) || (nPts <= 0)) {
            return null;
        }

        startLocs = new int[nPts];
        endLocs = ((tmpEnd == null) ? null : new int[nPts]);
        data = new float[nPts];
        features = (tmpFeatures == null ? null : new String[nPts]);

        try {
            System.arraycopy(tmpStart, startIdx, startLocs, 0, nPts);
            System.arraycopy(tmpData, startIdx, data, 0, nPts);
            if (endLocs != null) {
                System.arraycopy(tmpEnd, startIdx, endLocs, 0, nPts);
            }
            if (features != null) {
                System.arraycopy(tmpFeatures, startIdx, features, 0, nPts);
            }

        }
        catch (Exception exception) {
            log.error("Exception getting raw data tile", exception);
            return null;
        }


        return new DataTile(startLocs, endLocs, data, features);
    }


    public TrackType getTrackType() {
        try {
            return dataset.getType();
        }
        catch (Exception exception) {
            return TrackType.OTHER;
        }
    }


    @Override
    public boolean isLogNormalized() {
        return dataset.isLogNormalized();
    }


    public double getDataMax() {
        return dataset.getDataMax();
    }


    public double getDataMin() {
        return dataset.getDataMin();
    }

    @Override
    public int getLongestFeature(String chr) {
        return dataset.getLongestFeature(chr);
    }
}
