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

package org.broad.igv.feature.gff;

import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import htsjdk.tribble.Feature;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.track.FeatureSource;

import java.io.IOException;
import java.util.*;

/**
 * User: jacob
 * Date: 2012-Jun-22
 */
public class GFFFeatureSource implements org.broad.igv.track.FeatureSource {

    private static Logger log = Logger.getLogger(GFFFeatureSource.class);
    private  GFFCodec.Version gffVersion;

    private FeatureSource wrappedSource;

    public static boolean isGFF(String path) {
        String lowpath = path.toLowerCase();
        if (lowpath.endsWith(".gz")) {
            int idx = lowpath.length() - 3;
            lowpath = lowpath.substring(0, idx);
        }
        if (lowpath.endsWith(".txt")) {
            int idx = lowpath.length() - 4;
            lowpath = lowpath.substring(0, idx);
        }
        return lowpath.endsWith("gff3") || lowpath.endsWith("gvf") || lowpath.endsWith("gff") || lowpath.endsWith("gtf");
    }

    public static GFFCombiner getCombiner(GFFCodec.Version version) {
        return version == GFFCodec.Version.GFF3 ? new GFF3Combiner() : new GFF2Combiner();
    }

    public GFFFeatureSource(FeatureSource wrappedSource, GFFCodec.Version gffVersion) throws IOException {
        this.gffVersion = gffVersion;
        this.wrappedSource = wrappedSource;
    }

    @Override
    public void dispose() {
        wrappedSource.dispose();
    }

    @Override
    public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {

        // Expand start/end to be sure we get adjacent CDS/Exon features
        int expandedStart = Math.max(0, start - 2000000);
        long longEnd = (long) end + 2000000;   // Protect against overflow
        int expandedEnd = (int) Math.min(Integer.MAX_VALUE, longEnd);

        Iterator<Feature> rawIter = wrappedSource.getFeatures(chr, expandedStart, expandedEnd);
        GFFCombiner combiner = (getCombiner(gffVersion)).addFeatures(rawIter);

        // Now trim features not requested
        List<Feature> requestedFeatures = new ArrayList<>();
        for(Feature f : combiner.combineFeatures()) {
            if(f.getEnd() < start) continue;
            if(f.getStart() > end) break;
            requestedFeatures.add(f);
        }

        return new WrappedIterator(requestedFeatures.iterator());
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return wrappedSource.getCoverageScores(chr, start, end, zoom);
    }

    @Override
    public int getFeatureWindowSize() {
        return wrappedSource.getFeatureWindowSize();
    }

    @Override
    public void setFeatureWindowSize(int size) {
        wrappedSource.setFeatureWindowSize(size);
    }

}
