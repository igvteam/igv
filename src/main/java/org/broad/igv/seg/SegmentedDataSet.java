/*
 * sample       chrom   loc.start       loc.end num.mark        num.informative seg.mean
TCGA-02-0001-01C-01D-0183-04    1       554268  74674720        6892    6721    0.2077
TCGA-02-0001-01C-01D-0183-04    1       74693652        75110251        37      37      -0.2659
 */
package org.broad.igv.seg;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.util.*;


/**
 * @author jrobinso
 */
public class SegmentedDataSet {

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

    public SegmentedDataSet(Genome genome) {
        this.genome = genome;
    }

    public SegmentedDataSet(ResourceLocator locator, Genome genome) {
        //parser = locator.getPath().toLowerCase().endsWith(".gbench") ?
        //        new GBenchFileParser(locator) :
        //        new SegmentFileParser(locator);
        this.genome = genome;
        sortLists();


    }
    /**
     * Sort all segment lists.  This is done automatically after loading.
     */
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

        String chr = genome == null ? c : genome.getCanonicalChrName(c);

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
     * @param sample
     * @param chr
     * @return
     */
    public List<LocusScore> getSegments(String sample, String chr) {

        if(Globals.CHR_ALL.equals(chr)) {
            return getWholeGenomeScores(sample);
        }

        Map<String, List<LocusScore>> chrSegments = segments.get(sample);
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
    public double getDataMax() {
        return dataMax;
    }

    /**
     * Method description
     *
     * @param chr
     * @return
     */
    public double getDataMin() {
        return dataMin;
    }

    /**
     * Method description
     *
     * @param sample
     * @return
     */
    public List<LocusScore> getWholeGenomeScores(String sample) {
        List<LocusScore> wholeGenomeScores = wholeGenomeScoresCache.get(sample);
        if ((wholeGenomeScores == null) || wholeGenomeScores.isEmpty()) {
            wholeGenomeScores = new ArrayList(1000);
            for (String chr : genome.getLongChromosomeNames()) {
                List<LocusScore> chrSegments = getSegments(sample, chr);
                if (chrSegments != null) {
                    int lastgEnd = -1;
                    for (LocusScore score : chrSegments) {
                        Segment seg = (Segment) score;
                        int gStart = genome.getGenomeCoordinate(chr, seg.getStart());
                        int gEnd = genome.getGenomeCoordinate(chr, seg.getEnd());
                        if (gEnd > lastgEnd) {
                            wholeGenomeScores.add(new Segment(gStart, gStart, gEnd,
                                    gEnd, seg.getScore(), seg.getDescription()));
                        }
                        lastgEnd = gEnd;
                    }
                }
            }
            wholeGenomeScoresCache.put(sample, wholeGenomeScores);
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

    public LocusScore getSegmentAt(String sample, String chr, double position, ReferenceFrame frame) {

        List<LocusScore> scores = getSegments(sample, chr);

        if (scores == null) {
            return null;
        } else {
            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            int buffer = (int) (2 * bpPerPixel);    /* * */
            return (LocusScore) FeatureUtils.getFeatureAt(position, buffer, scores);
        }
    }


}
