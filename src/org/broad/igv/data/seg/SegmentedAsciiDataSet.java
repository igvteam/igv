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

    //SegFileParser parser;
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
     *
     */
    public void addSegment(String heading, String c, int start, int end, float value, String desc) {

        String chr = genome == null ? c : genome.getChromosomeAlias(c);

        Map<String, List<LocusScore>> chrSegments = segments.get(heading);
        if (chrSegments == null) {
            headings.add(heading);
            chrSegments = new HashMap();
            segments.put(heading, chrSegments);
        }

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

            // Compute the smallest concievable feature that could be viewed on the
            // largest screen.  Be conservative.   The smallest feature is one at
            // the screen resolution scale in <chr units> / <pixel>
            double minFeatureSize = 0; // ((double) genome.getLength()) / (maxScreenSize * locationUnit);

            long offset = 0;
            wholeGenomeScores = new ArrayList(1000);
            for (String chr : genome.getLongChromosomeNames()) {
                List<LocusScore> chrSegments = getSegments(heading, chr);
                if (chrSegments != null) {
                    int lastgEnd = -1;
                    for (LocusScore score : chrSegments) {
                        Segment seg = (Segment) score;
                        int gStart = genome.getGenomeCoordinate(chr, seg.getStart());
                        int gEnd = genome.getGenomeCoordinate(chr, seg.getEnd());
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
        if (trackProperties == null) {
            trackProperties = new TrackProperties();
        }
        return trackProperties;
    }


}
