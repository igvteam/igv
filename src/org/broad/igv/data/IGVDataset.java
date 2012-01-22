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

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ObjectCache;
import org.broad.igv.util.ResourceLocator;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class IGVDataset implements Dataset {

    private String name;

    private TrackType type = TrackType.OTHER;

    private boolean logNormalized;
    private String[] dataHeadings;
    private Map<String, ChromosomeSummary> chromosomeSummaries = new LinkedHashMap();
    private GenomeSummaryData genomeSummary;
    private IGVDatasetParser parser;
    private ObjectCache<String, ChromosomeData> chromsomeDataCache = new ObjectCache(30);
    private float dataMin;
    private float dataMax;
    TrackProperties trackProperties = new TrackProperties();
    private Map<String, Integer> longestFeatureMap;

    public IGVDataset(ResourceLocator locator, Genome genome, IGV igv) {

        parser = new IGVDatasetParser(locator, genome, igv);

        List<ChromosomeSummary> summaries = parser.scan(this);

        if (summaries == null || summaries.size() == 0)
            throw new RuntimeException("Could not find any chromosomes in the dataset on the genome(" + genome.getId() + ")");

        for (ChromosomeSummary summary : summaries) {
            chromosomeSummaries.put(summary.getName(), summary);
        }
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void setGenomeSummary(GenomeSummaryData genomeSummary) {
        this.genomeSummary = genomeSummary;
    }

    public void setTrackType(TrackType type) {
        this.type = type;
    }

    public TrackType getType() {
        return type;
    }

    public String[] getChromosomes() {
        return chromosomeSummaries.keySet().toArray(new String[]{});
    }

    public void setDataHeadings(String[] dataHeadings) {
        this.dataHeadings = dataHeadings;
    }

    public String[] getTrackNames() {
        return dataHeadings;
    }

    public int[] getStartLocations(String chr) {
        ChromosomeData cd = getChromosomeData(chr);
        return cd == null ? null : cd.getStartLocations();
    }

    public int[] getEndLocations(String chr) {
        ChromosomeData cd = getChromosomeData(chr);
        return cd == null ? null : cd.getEndLocations();
    }

    public String[] getFeatureNames(String chr) {
        ChromosomeData cd = getChromosomeData(chr);
        return cd == null ? null : cd.getProbes();
    }

    public float[] getData(String heading, String chr) {
        ChromosomeData cd = getChromosomeData(chr);
        return cd == null ? null : cd.getData(heading);
    }

    /**
     * Get the data for all samples (tracks) for the given chromosome.
     * <p/>
     * This method is synchronized to insure that the data for a chromosome
     * is only loaded once.
     *
     * @param chr
     * @return
     */
    private synchronized ChromosomeData getChromosomeData(String chr) {
        ChromosomeData cd = chromsomeDataCache.get(chr);
        if (cd == null) {
            ChromosomeSummary sum = chromosomeSummaries.get(chr);
            if (sum == null) {
                //todo -- throw exception
                return null;
            }
            //synchronized (sum) {

            cd = parser.loadChromosomeData(sum, dataHeadings);
            chromsomeDataCache.put(chr, cd);
            //}

        }
        return cd;
    }

    public GenomeSummaryData getGenomeSummary() {
        return genomeSummary;
    }

    public void setLogNormalized(boolean logNormalized) {
        this.logNormalized = logNormalized;
    }

    public boolean isLogNormalized() {
        return logNormalized;
    }

    public float getDataMin() {
        return dataMin;
    }

    public float getDataMax() {
        return dataMax;
    }

    public void setDataMin(float dataMin) {
        this.dataMin = dataMin;
    }

    public void setDataMax(float dataMax) {
        this.dataMax = dataMax;
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }


    public Integer getLongestFeature(String chr) {
        return longestFeatureMap == null ? 1000 :
                longestFeatureMap.containsKey(chr) ? longestFeatureMap.get(chr) : 1;
    }

    public void setLongestFeatureMap(Map<String, Integer> longestFeatureMap) {
        this.longestFeatureMap = longestFeatureMap;
    }
}
