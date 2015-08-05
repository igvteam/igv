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

package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.EncodePeakFeature;
import org.broad.igv.feature.FeatureUtils;
import htsjdk.tribble.Feature;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 1/14/14
 *         Time: 2:04 PM
 */
public class CursorRegion implements Feature {

    String chr;  //
    int location;  // Genomic location in BP.  Zero based coords

    public CursorRegion(String chr, int location) {
        this.chr = chr;
        this.location = location;
    }

    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public int getStart() {
        return location;
    }

    @Override
    public int getEnd() {
        return location + 1;
    }

    public int getLocation() {
        return location;
    }

    public double getScore(CursorTrack t,  int frameBPWidth) {

        int longest = t.getLongestFeatureLength(chr);
        List<BasicFeature> features = t.getFeatures(chr);

        double score = -1;
        if (features == null) return score;

        FeatureIterator iter = getFeatureIterator(features, longest, frameBPWidth);
        while (iter.hasNext()) {
            BasicFeature f = iter.next();
            score = Math.max(t.getSignal(f), score);
        }
        return score;
    }

    public FeatureIterator getFeatureIterator(List<BasicFeature> features, int longest, int frameBPWidth) {
        return new FeatureIterator(features, longest, frameBPWidth);
    }

    class FeatureIterator implements Iterator<BasicFeature> {

        int bpStart;
        int bpEnd;
        int longest;
        BasicFeature nextFeature;
        Iterator<BasicFeature> features;

        FeatureIterator(List<BasicFeature> features, int longest, int frameBPWidth) {
            this.longest = longest;
            bpStart = location - frameBPWidth / 2;
            bpEnd = location + frameBPWidth/2;

            int s0 = longest < 0 ? 0 : bpStart - longest;
            int startIdx = FeatureUtils.getIndexBefore(s0, features);
            this.features = (features.subList(startIdx, features.size())).iterator();
            advance();


        }

        private void advance() {
            // We have to look back at least as far as the longest feature for overlap.
            nextFeature = null;
            while (features.hasNext()) {
                BasicFeature f = features.next();
                if (f.getStart() >= bpEnd) {
                    break;
                }

                if (f.getEnd() >= bpStart && f.getStart() < bpEnd) {
                    nextFeature = f;
                    break;
                }
            }
        }

        @Override
        public boolean hasNext() {
            return nextFeature != null;
        }

        @Override
        public BasicFeature next() {
            BasicFeature retValue = nextFeature;
            advance();
            return retValue;
        }

        @Override
        public void remove() {

        }
    }

}
