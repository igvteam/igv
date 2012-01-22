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

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;

import java.util.HashMap;
import java.util.List;
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
    public DatasetDataSource(String trackId, Dataset dataset, Genome genome) {
        super(genome);
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
    protected DataTile getRawData(String chr, int startLocation, int endLocation) {

        if (chr.equals(Globals.CHR_ALL) && genomeSummaryData != null) {
            int[] startLocs = genomeSummaryData.getLocations();
            int[] endLocs = null;
            float[] data = genomeSummaryData.getData(trackId);
            return new DataTile(startLocs, endLocs, data, null);
        } else {
            int[] startLocs = dataset.getStartLocations(chr);;
            int[] endLocs = dataset.getEndLocations(chr);
            float[] data =  dataset.getData(trackId, chr);
            String[] features = dataset.getFeatureNames(chr);

            if (startLocs == null|| (data == null) || data.length == 0) {
                return null;
            }

            return new DataTile(startLocs, endLocs, data, features);
        }
    }

    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
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
