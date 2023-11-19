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

package org.broad.igv.ucsc.bb;

import htsjdk.tribble.Feature;
import org.apache.commons.math3.stat.StatUtils;
import org.broad.igv.ucsc.bb.codecs.BBCodec;
import org.broad.igv.ucsc.bb.codecs.BBCodecFactory;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ucsc.UnsignedByteBuffer;

import java.io.IOException;
import java.util.*;

/**
 * A hybrid source, implements both DataSource and FeatureSource.
 *
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BBDataSource implements FeatureSource, DataSource {

    final int screenWidth = 1000; // TODO use actual screen width


    Collection<WindowFunction> availableWindowFunctions =
            Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max, WindowFunction.none);

    BBFile reader;

    // Feature visibility window (for bigBed)
    private int featureVisiblityWindow = -1;

    private Map<WindowFunction, List<LocusScore>> wholeGenomeScores;

    // Lookup table to support chromosome aliasing.
    private Map<String, String> chrNameMap = new HashMap();
    private double dataMin = 0;
    private double dataMax = 100;

    BBCodec bedCodec;

    public BBDataSource(BBFile reader, Genome genome) throws IOException {
        super();

        this.reader = reader;

        this.wholeGenomeScores = new HashMap<>();

        if (reader.type == BBFile.Type.BIGWIG) initMinMax();

        if (reader.type == BBFile.Type.BIGBED) {
            String autosql = reader.autosql;
            int definedFieldCount = reader.header.definedFieldCount;
            bedCodec = BBCodecFactory.getCodec(autosql, definedFieldCount);
        }


        // Assume 1000 pixel screen, pick visibility level to be @ highest resolution zoom.
        // NOTE: this is only used by feature tracks (bigbed sources)
        // TODO -- something smarter, like scaling by actual density
        //if (levels != null && levels.getZoomHeaderCount() > 0) {
        //    BBZoomLevelHeader firstLevel = levels.getZoomLevelHeaders().get(0); // Highest res
        featureVisiblityWindow = 0;
        //}

//        if (genome != null) {
//            Collection<String> chrNames = reader.getChromosomeNames();
//            for (String chr : chrNames) {
//                String igvChr = genome.getCanonicalChrName(chr);
//                if (igvChr != null && !igvChr.equals(chr)) {
//                    chrNameMap.put(igvChr, chr);
//                }
//            }
//        }

    }

    @Override
    public int getFeatureWindowSize() {
        return featureVisiblityWindow;
    }

    /**
     * Set the "min" and "max" from 1MB resolutiond data.  Read a maximum of 10,000 points for this
     */
    private void initMinMax() {

        final int oneMB = 1000000;
        //     final BBZoomLevelHeader zoomLevelHeader = getZoomLevelForScale(oneMB);

        int nValues = 0;
        double[] values = new double[10000];

//        if (zoomLevelHeader == null) {
//            List<String> chrNames = reader.getChromosomeNames();
//            for (String chr : chrNames) {
//                BigWigIterator iter = reader.getBigWigIterator(chr, 0, chr, Integer.MAX_VALUE, false);
//                while (iter.hasNext()) {
//                    WigItem item = iter.next();
//                    values[nValues++] = item.getWigValue();
//                    if (nValues >= 10000) break;
//                }
//            }
//        } else {
//
//            int z = zoomLevelHeader.getZoomLevel();
//            ZoomLevelIterator zlIter = reader.getZoomLevelIterator(z);
//            if (zlIter.hasNext()) {
//                while (zlIter.hasNext()) {
//                    ZoomDataRecord rec = zlIter.next();
//                    values[nValues++] = (rec.getMeanVal());
//                    if (nValues >= 10000) {
//                        break;
//                    }
//                }
//            }
//        }

        if (nValues > 0) {
            dataMin = StatUtils.percentile(values, 0, nValues, 10);
            dataMax = StatUtils.percentile(values, 0, nValues, 90);
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

    @Override
    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {
        return null;
    }

    public TrackType getTrackType() {
        return TrackType.OTHER;
    }

    @Override
    public void setWindowFunction(WindowFunction statType) {

    }

    public boolean isLogNormalized() {
        return false;
    }

    @Override
    public WindowFunction getWindowFunction() {
        return null;
    }


    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return availableWindowFunctions;
    }

    @Override
    public void dispose() {

    }

    @Override
    public void close() {
        // super.dispose();
        if (reader != null) {
            //     reader.close();
        }
    }

    // Feature interface follows ------------------------------------------------------------------------

    public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {

        List<byte[]> chunks = this.reader.getLeafChunks(chr, start, chr, end, 0);

        List features = new ArrayList<>();
        for (byte[] c : chunks) {
            UnsignedByteBuffer bb = UnsignedByteBuffer.wrap(c, reader.byteOrder);
            while (bb.remaining() > 0) {

                int chromId = bb.getInt();
                int chromStart = bb.getInt();
                int chromEnd = bb.getInt();
                String restOfFields = bb.getString();
                //String chr = reader.getChrForId(chromId);  // Might differ due to aliasing

                final BedData bedData = new BedData(chr, chromStart, chromEnd, restOfFields);
                final BasicFeature feature = bedCodec.decode(bedData);
                features.add(feature);
            }
        }

        return new FeatureIterator(features, start, end);
    }

    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    static class FeatureIterator implements Iterator<Feature> {

        List<Feature> features;
        int idx;
        int start;
        int end;

        Feature next;

        public FeatureIterator(List<Feature> features, int start, int end) {
            this.features = features;
            this.start = start;
            this.end = end;
            advance();
        }

        private void advance() {
            while(idx < features.size()) {
                next = features.get(idx++);
                if(next.getStart() > end) {
                    next = null;
                    break;   // Done
                } else if(next.getEnd() >= start) {
                    break;
                }
            }
        }

        @Override
        public boolean hasNext() {
            return next != null;
        }

        @Override
        public Feature next() {
            Feature retValue = next;
            advance();
            return retValue;
        }
    }


}
