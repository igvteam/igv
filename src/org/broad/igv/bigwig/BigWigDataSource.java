/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.bigwig;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.Globals;
import org.broad.igv.bbfile.*;
import org.broad.igv.data.AbstractDataSource;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataTile;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.util.collections.IntArrayList;
import org.broad.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * A hybrid source, implements both DataSource and FeatureSource.   Way of the future?
 *
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BigWigDataSource extends AbstractDataSource implements FeatureSource {

    final int screenWidth = 1000; // TODO use actual screen width


    Collection<WindowFunction> availableWindowFunctions =
            Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max);
    WindowFunction windowFunction = WindowFunction.mean;

    BBFileReader reader;
    private BBZoomLevels levels;

    // Feature visibility window (for bigBed)
    private int featureVisiblityWindow = -1;

    private List<LocusScore> wholeGenomeScores;

    // Lookup table to support chromosome aliasing.
    private Map<String, String> chrNameMap = new HashMap();

    private RawDataInterval currentInterval = null;


    private double dataMin = 0;
    private double dataMax = 100;

    public BigWigDataSource(BBFileReader reader, Genome genome) throws IOException {
        super(genome);

        this.reader = reader;
        levels = reader.getZoomLevels();
        initMinMax();

        // Assume 1000 pixel screen, pick visibility level to be @ highest resolution zoom.
        // TODO -- something smarter, like scaling by actual density
        if (levels.getZoomHeaderCount() > 0) {
            BBZoomLevelHeader firstLevel = levels.getZoomLevelHeaders().get(0); // Highest res
            featureVisiblityWindow = firstLevel.getReductionLevel() * 2000;
        }

        if (genome != null) {
            Collection<String> chrNames = reader.getChromosomeNames();
            for (String chr : chrNames) {
                String igvChr = genome.getChromosomeAlias(chr);
                if (igvChr != null && !igvChr.equals(chr)) {
                    chrNameMap.put(igvChr, chr);
                }
            }
        }
    }


    /**
     * Set the "min" and "max" from 1MB resolutiond data.  Read a maximum of 10,000 points for this
     */
    private void initMinMax() {
        final int oneMB = 1000000;
        int z = getZoomLevelForScale(oneMB).getZoomLevel();

        ZoomLevelIterator zlIter = reader.getZoomLevelIterator(z);
        int n = 0;
        double[] values = new double[10000];
        if (zlIter.hasNext()) {
            while (zlIter.hasNext()) {
                ZoomDataRecord rec = zlIter.next();
                values[n] = (rec.getMeanVal());
                n++;
                if (n >= 10000) {
                    break;
                }
            }
            dataMin = StatUtils.percentile(values, 0, n, 10);

            // Peg the scale at 100, this seems arbitrary but it maintains some compatibility with earlier releases
            // for 'large-valued' data, but corrects scaling problems for small-valued data.
            dataMax = Math.min(100, StatUtils.percentile(values, 0, n, 90));
        } else {
            dataMin = 0;
            dataMax = 100;
        }
    }

    public double getDataMax() {
        return dataMax;
    }

    public double getDataMin() {
        return dataMin;
    }

    public TrackType getTrackType() {
        return TrackType.OTHER;
    }

    public void setWindowFunction(WindowFunction statType) {
        // Invalidate caches
        wholeGenomeScores = null;

        this.windowFunction = statType;
    }

    public boolean isLogNormalized() {
        return false;
    }

    public void refreshData(long timestamp) {

    }

    @Override
    public int getLongestFeature(String chr) {
        return 0;
    }

    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return availableWindowFunctions;
    }

    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int start, int end, int zoom) {

        if (chr.equals(Globals.CHR_ALL)) {
            return getWholeGenomeScores();
        } else {
            return getZoomSummaryScores(chr, start, end, zoom);
        }
    }


    /**
     * Return the zoom level that most closely matches the given resolution.  Resolution is in BP / Pixel.
     *
     * @param resolution
     * @return
     */
    private BBZoomLevelHeader getZoomLevelForScale(double resolution) {
        final ArrayList<BBZoomLevelHeader> headers = levels.getZoomLevelHeaders();

        BBZoomLevelHeader lastLevel = null;
        for (BBZoomLevelHeader zlHeader : headers) {
            int reductionLevel = zlHeader.getReductionLevel();
            if (reductionLevel > resolution) {
                return lastLevel == null ? zlHeader : lastLevel;
            }
            lastLevel = zlHeader;
        }
        return headers.get(headers.size() - 1);

    }

    private BBZoomLevelHeader getLowestResolutionLevel() {
        final ArrayList<BBZoomLevelHeader> headers = levels.getZoomLevelHeaders();
        return headers.get(headers.size() - 1);
    }

    protected List<LocusScore> getZoomSummaryScores(String chr, int start, int end, int zoom) {

        Chromosome c = genome.getChromosome(chr);
        if (c == null) return null;

        double nBins = Math.pow(2, zoom);

        double scale = c.getLength() / (nBins * 700);

        BBZoomLevelHeader zlHeader = getZoomLevelForScale(scale);
        int bbLevel = zlHeader.getZoomLevel();
        int reductionLevel = zlHeader.getReductionLevel();


        // If we are at the highest precomputed resolution compare to the requested resolution.  If they differ
        // by more than a factor of 2 compute "on the fly"
        String tmp = chrNameMap.get(chr);
        String querySeq = tmp == null ? chr : tmp;

        if (reader.isBigBedFile() || bbLevel > 1 || (bbLevel == 1 && (reductionLevel / scale) < 2)) {
            ArrayList<LocusScore> scores = new ArrayList(1000);
            ZoomLevelIterator zlIter = reader.getZoomLevelIterator(bbLevel, querySeq, start, querySeq, end, false);
            while (zlIter.hasNext()) {
                ZoomDataRecord rec = zlIter.next();

                float v = getValue(rec);
                BasicScore bs = new BasicScore(rec.getChromStart(), rec.getChromEnd(), v);
                scores.add(bs);
            }
            return scores;

        } else {
            // No precomputed scores for this resolution level
            return null;
        }
    }

    private float getValue(ZoomDataRecord rec) {
        float v;
        switch (windowFunction) {
            case min:
                v = rec.getMinVal();
                break;
            case max:
                v = rec.getMaxVal();
                break;
            default:
                v = rec.getMeanVal();

        }
        return v;
    }


    @Override
    protected synchronized DataTile getRawData(String chr, int start, int end) {

        if (chr.equals(Globals.CHR_ALL)) {
            return null;
        }


        if (currentInterval != null && currentInterval.contains(chr, start, end)) {
            return currentInterval.tile;
        }

        // TODO -- fetch data directly in arrays to avoid creation of multiple "WigItem" objects?
        IntArrayList startsList = new IntArrayList(100000);
        IntArrayList endsList = new IntArrayList(100000);
        FloatArrayList valuesList = new FloatArrayList(100000);

        String chrAlias = chrNameMap.containsKey(chr) ? chrNameMap.get(chr) : chr;
        Iterator<WigItem> iter = reader.getBigWigIterator(chrAlias, start, chrAlias, end, false);

        while (iter.hasNext()) {
            WigItem wi = iter.next();
            startsList.add(wi.getStartBase());
            endsList.add(wi.getEndBase());
            valuesList.add(wi.getWigValue());
        }

        DataTile tile = new DataTile(startsList.toArray(), endsList.toArray(), valuesList.toArray(), null);
        currentInterval = new RawDataInterval(chr, start, end, tile);

        return tile;

    }


    private List<LocusScore> getWholeGenomeScores() {

        if (genome.getHomeChromosome().equals(Globals.CHR_ALL)) {
            if (wholeGenomeScores == null) {
                double scale = genome.getNominalLength() / screenWidth;
                wholeGenomeScores = new ArrayList<LocusScore>();

                for (String chrName : genome.getLongChromosomeNames()) {
                    Chromosome chr = genome.getChromosome(chrName);

                    BBZoomLevelHeader lowestResHeader = this.getZoomLevelForScale(scale);

                    int lastGenomeEnd = -1;
                    int end = chr.getLength();

                    String tmp = chrNameMap.get(chrName);
                    String querySeq = tmp == null ? chrName : tmp;

                    ZoomLevelIterator zlIter = reader.getZoomLevelIterator(
                            lowestResHeader.getZoomLevel(), querySeq, 0, querySeq, end, false);

                    while (zlIter.hasNext()) {
                        ZoomDataRecord rec = zlIter.next();

                        float value = getValue(rec);
                        int genomeStart = genome.getGenomeCoordinate(chrName, rec.getChromStart());
                        if (genomeStart < lastGenomeEnd) {
                            continue;
                        }

                        int genomeEnd = genome.getGenomeCoordinate(chrName, rec.getChromEnd());
                        wholeGenomeScores.add(new BasicScore(genomeStart, genomeEnd, value));
                        lastGenomeEnd = genomeEnd;
                    }
                }


            }
            return wholeGenomeScores;
        } else {
            return null;
        }

    }


    // Feature interface follows ------------------------------------------------------------------------

    public Iterator getFeatures(String chr, int start, int end) throws IOException {

        String tmp = chrNameMap.get(chr);
        String querySeq = tmp == null ? chr : tmp;
        BigBedIterator bedIterator = reader.getBigBedIterator(querySeq, start, chr, end, false);
        return new WrappedIterator(bedIterator);
    }

    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        String tmp = chrNameMap.get(chr);
        String querySeq = tmp == null ? chr : tmp;
        return this.getSummaryScoresForRange(querySeq, start, end, zoom);
    }

    public int getFeatureWindowSize() {
        return this.featureVisiblityWindow;
    }

    public void setFeatureWindowSize(int size) {
        this.featureVisiblityWindow = size;
    }

    public Class getFeatureClass() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public static class WrappedIterator implements Iterator<Feature> {

        BigBedIterator bedIterator;

        public WrappedIterator(BigBedIterator bedIterator) {
            this.bedIterator = bedIterator;
        }

        public boolean hasNext() {
            return bedIterator.hasNext();  //To change body of implemented methods use File | Settings | File Templates.
        }

        public Feature next() {
            BedFeature feat = bedIterator.next();
            BasicFeature feature = new BasicFeature(feat.getChromosome(), feat.getStartBase(), feat.getEndBase());
            String[] restOfFields = feat.getRestOfFields();
            if (restOfFields != null && restOfFields.length > 0) {
                decode(feature, restOfFields);
            }
            return feature;

        }

        public void remove() {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }


    //  End FeatureSource interface ----------------------------------------------------------------------

    static class RawDataInterval {
        String chr;
        int start;
        int end;
        DataTile tile;

        RawDataInterval(String chr, int start, int end, DataTile tile) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.tile = tile;
        }

        public boolean contains(String chr, int start, int end) {
            return chr.equals(this.chr) && start >= this.start && end <= this.end;
        }
    }


    /////////// Decoder for BED features


    private static void decode(BasicFeature feature, String[] restOfFields) {

        int tokenCount = restOfFields.length;


        // The rest of the columns are optional.  Stop parsing upon encountering
        // a non-expected value

        // Name
        if (tokenCount > 0) {
            String name = restOfFields[0].replaceAll("\"", "");
            feature.setName(name);
            feature.setIdentifier(name);
        }

        // Score
        if (tokenCount > 1) {
            try {
                float score = Float.parseFloat(restOfFields[1]);
                feature.setScore(score);
            } catch (NumberFormatException numberFormatException) {

                // Unexpected, but does not invalidate the previous values.
                // Stop parsing the line here but keep the feature
                // Don't log, would just slow parsing down.
                return;
            }
        }

        // Strand
        if (tokenCount > 2) {
            String strandString = restOfFields[2].trim();
            char strand = (strandString.length() == 0)
                    ? ' ' : strandString.charAt(0);

            if (strand == '-') {
                feature.setStrand(Strand.NEGATIVE);
            } else if (strand == '+') {
                feature.setStrand(Strand.POSITIVE);
            } else {
                feature.setStrand(Strand.NONE);
            }
        }

        if (tokenCount > 5) {
            String colorString = restOfFields[5];
            feature.setColor(ColorUtilities.stringToColor(colorString));
        }

        // Coding information is optional
        if (tokenCount > 8) {
            Strand strand = feature.getStrand();
            int cdStart = Integer.parseInt(restOfFields[3]);
            int cdEnd = Integer.parseInt(restOfFields[4]);

            int exonCount = Integer.parseInt(restOfFields[6]);
            String[] exonSizes = Globals.commaPattern.split(restOfFields[7]);
            String[] startsBuffer = Globals.commaPattern.split(restOfFields[8]);

            int exonNumber = (strand == Strand.NEGATIVE ? exonCount : 1);

            int start = feature.getStart();
            String chr = feature.getChr();
            if (startsBuffer.length == exonSizes.length) {
                for (int i = 0; i < startsBuffer.length; i++) {
                    int exonStart = start + Integer.parseInt(startsBuffer[i]);
                    int exonEnd = exonStart + Integer.parseInt(exonSizes[i]);
                    Exon exon = new Exon(chr, exonStart, exonEnd, strand);
                    exon.setCodingStart(cdStart);
                    exon.setCodingEnd(cdEnd);
                    exon.setNumber(exonNumber);
                    feature.addExon(exon);

                    if (strand == Strand.NEGATIVE) {
                        exonNumber--;
                    } else {
                        exonNumber++;
                    }
                }
            }
        }
    }


}
