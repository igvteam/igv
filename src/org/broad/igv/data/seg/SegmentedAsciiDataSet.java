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
 * sample       chrom   loc.start       loc.end num.mark        num.informative seg.mean
TCGA-02-0001-01C-01D-0183-04    1       554268  74674720        6892    6721    0.2077
TCGA-02-0001-01C-01D-0183-04    1       74693652        75110251        37      37      -0.2659
 */
package org.broad.igv.data.seg;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ResourceLocator;

import java.util.*;


/**
 * @author jrobinso
 */
public class SegmentedAsciiDataSet implements SegmentedDataSet {

    SegFileParser parser;
    TrackType trackType = TrackType.COPY_NUMBER;
    float dataMax = -Float.MAX_VALUE;
    float dataMin = Float.MAX_VALUE;
    /**
     * Assume data is non-log value until suggested otherwise by the precense
     * of negative numbers.  TODO This is a fragile assumption, the user should
     * input this information directly.
     */
    private boolean logNormalized = false;
    /**
     * Map of [heading ->  [chr -> [list of chrSegments]]]
     */
    private Map<String, Map<String, List<LocusScore>>> segments = new HashMap();
    /**
     * Set of chromosomes represented in this dataset
     */
    private Set<String> chromosomes = new HashSet();
    private List<String> headings = new ArrayList();
    private Map<String, List<LocusScore>> wholeGenomeScoresCache = new HashMap();
    private long lastRefreshTime = 0;
    private TrackProperties trackProperties;
    Genome genome;

    public SegmentedAsciiDataSet(Genome genome) {
        this.genome = genome;
    }

    public SegmentedAsciiDataSet(ResourceLocator locator, Genome genome) {
        //parser = locator.getPath().toLowerCase().endsWith(".gbench") ?
        //        new GBenchFileParser(locator) :
        //        new SegmentFileParser(locator);
        this.genome = genome;
        parser = new SegmentFileParser(locator);
        parser.loadSegments(this, genome);
        sortLists();


    }


    public void sortLists() {
        for (Map.Entry<String, Map<String, List<LocusScore>>> sampleEntry : segments.entrySet()) {
            for (Map.Entry<String, List<LocusScore>> chrEntry : sampleEntry.getValue().entrySet()) {
                List<LocusScore> tmp = chrEntry.getValue();
                FeatureUtils.sortFeatureList(tmp);
            }
        }
    }

    /**
     * Method description
     *
     * @param timestamp
     */
    public synchronized void refreshData(long timestamp) {
        if (timestamp > lastRefreshTime) {
            dataMax = -Float.MAX_VALUE;
            dataMin = Float.MAX_VALUE;
            logNormalized = false;
            segments = new HashMap();
            headings = new ArrayList();
            wholeGenomeScoresCache.clear();
            parser.loadSegments(this, genome);
            lastRefreshTime = timestamp;
        }

    }

    /**
     *
     */
    public void addSegment(String heading, String c, int start, int end, float value, String desc) {

        Map<String, List<LocusScore>> chrSegments = segments.get(heading);
        if (chrSegments == null) {
            headings.add(heading);
            chrSegments = new HashMap();
            segments.put(heading, chrSegments);
        }

        String chr = genome == null ? c : genome.getChromosomeAlias(c);
        List<LocusScore> segmentList = chrSegments.get(chr);
        if (segmentList == null) {
            segmentList = new ArrayList<LocusScore>();
            chrSegments.put(chr, segmentList);
        }
        segmentList.add(new Segment(start, start, end, end, value, desc));
        dataMax = Math.max(dataMax, value);
        dataMin = Math.min(dataMin, value);
        if (value < 0) {
            logNormalized = true;
        }

        chromosomes.add(chr);

    }


    /**
     * Method description
     *
     * @return
     */
    public Set<String> getChromosomes() {
        return chromosomes;
    }

    /**
     * Method description
     *
     * @param heading
     * @param chr
     * @return
     */
    public List<LocusScore> getSegments(String heading, String chr) {
        Map<String, List<LocusScore>> chrSegments = segments.get(heading);
        return (chrSegments == null) ? null : chrSegments.get(chr);
    }

    /**
     * Method description
     *
     * @return
     */
    public List<String> getDataHeadings() {
        return headings;
    }


    public List<String> getSampleNames() {
        return headings;
    }

    /**
     * Method description
     *
     * @return
     */
    public TrackType getType() {
        return trackType;
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isLogNormalized() {
        return logNormalized;
    }

    /**
     * Method description
     *
     * @param chr
     * @return
     */
    public double getDataMax(String chr) {
        return dataMax;
    }

    /**
     * Method description
     *
     * @param chr
     * @return
     */
    public double getDataMin(String chr) {
        return dataMin;
    }

    /**
     * Method description
     *
     * @param heading
     * @return
     */
    public List<LocusScore> getWholeGenomeScores(String heading) {


        List<LocusScore> wholeGenomeScores = wholeGenomeScoresCache.get(heading);
        if ((wholeGenomeScores == null) || wholeGenomeScores.isEmpty()) {
            int locationUnit = 1000;

            // Compute the smallest concievable feature that could be viewed on the
            // largest screen.  Be conservative.   The smallest feature is one at
            // the screen resolution scale in <chr units> / <pixel>
            int maxScreenSize = 4000;
            double minFeatureSize = 0; // ((double) genome.getLength()) / (maxScreenSize * locationUnit);

            long offset = 0;
            wholeGenomeScores = new ArrayList(1000);
            for (String chr : genome.getChromosomeNames()) {
                List<LocusScore> chrSegments = getSegments(heading, chr);
                if (chrSegments != null) {
                    int lastgEnd = -1;
                    for (LocusScore score : chrSegments) {
                        Segment seg = (Segment) score;
                        int gStart = (int) ((offset + score.getStart()) / locationUnit);
                        int gEnd = (int) ((offset + score.getEnd()) / locationUnit);
                        if ((gEnd - gStart) > minFeatureSize) {
                            wholeGenomeScores.add(new Segment(gStart, gStart, gEnd,
                                    gEnd, seg.getScore(), seg.getDescription()));
                        }
                    }

                }
                offset += genome.getChromosome(chr).getLength();
            }
            wholeGenomeScoresCache.put(heading, wholeGenomeScores);
        }
        return wholeGenomeScores;

    }

    /**
     * Method description
     *
     * @return
     */
    public TrackType getTrackType() {
        return trackType;
    }

    /**
     * Method description
     *
     * @param trackType
     */
    public void setTrackType(TrackType trackType) {
        this.trackType = trackType;
    }


    public void setTrackProperties(TrackProperties props) {
        this.trackProperties = props;
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }


}
