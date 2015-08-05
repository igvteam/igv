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
     *
     * @param trackId
     * @param dataset
     * @param genome
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
                        Map<String, float[]> dMap = new HashMap<String, float[]>();
                        dMap.put(trackId, dataset.getData(trackId, chr.getName()));
                        genomeSummaryData.addData(chr.getName(), startLocations, dMap);
                    }
                }
            }
        }
    }


    @Override
    protected DataTile getRawData(String chr, int startLocation, int endLocation) {

        if (chr.equals(Globals.CHR_ALL) && genomeSummaryData != null) {
            int[] startLocs = genomeSummaryData.getLocations();
            int[] endLocs = null;
            float[] data = genomeSummaryData.getData(trackId);
            return new DataTile(startLocs, endLocs, data, null);
        } else {
            int[] startLocs = dataset.getStartLocations(chr);
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
        return null;
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
