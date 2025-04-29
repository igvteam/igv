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

import htsjdk.tribble.NamedFeature;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.WindowFunction;

import java.io.IOException;
import java.util.*;

/**
 * A hybrid source, implements both DataSource and FeatureSource.
 *
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BBFeatureSource implements FeatureSource {

    private static Logger log = LogManager.getLogger(BBFeatureSource.class);
    private final Genome genome;

    Collection<WindowFunction> availableWindowFunctions =
            Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max, WindowFunction.none);

    BBFile reader;

    // Feature visibility window (for bigBed)
    private int featureVisiblityWindow = -1;

    private Map<WindowFunction, List<LocusScore>> wholeGenomeScores;

    public BBFeatureSource(BBFile reader, Genome genome) throws IOException {

        super();
        this.genome = genome;
        this.reader = reader;

        // viz window to load average of 10,000 features per screen
        featureVisiblityWindow = reader.getFeatureDensity() > 0 ? (int) (10000 / reader.getFeatureDensity()) : -1;

    }

    @Override
    public int getFeatureWindowSize() {
        return featureVisiblityWindow;
    }

    @Override
    public void close() {
        // super.dispose();
        if (reader != null) {
            //     reader.close();
        }
    }

    /**
     * Return an iterator for features spanning the genomic range.   If no features span the range returns an empy list.
     *
     * @param chr
     * @param start
     * @param end
     * @return
     * @throws IOException
     */
    public Iterator<BasicFeature> getFeatures(String chr, int start, int end) throws IOException {

        long rTreeOffset = reader.getHeader().fullIndexOffset;
        Integer chrIdx = reader.getIdForChr(chr);
        if(chrIdx == null) {
            return Collections.emptyIterator();
        } else {
            List<byte[]> chunks = this.reader.getLeafChunks(chrIdx, start, chrIdx, end, rTreeOffset);
            List features = new ArrayList<>();
            for (byte[] chunk : chunks) {
                features.addAll(reader.decodeFeatures(chr, chunk, chrIdx, start, end));
            }
            return new FeatureIterator(features, start, end);
        }
    }

    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    public boolean isSearchable() {
        return reader.isSearchable();
    }

    @Override
    public List<BasicFeature> search(String term) {
        try {
            return reader.search(term);
        } catch (IOException e) {
            log.error("Error searching for: " + term, e);
            return null;
        }
    }

    static class FeatureIterator implements Iterator<BasicFeature> {

        List<BasicFeature> features;
        int idx;
        int start;
        int end;

        BasicFeature next;

        public FeatureIterator(List<BasicFeature> features, int start, int end) {
            this.features = features;
            this.start = start;
            this.end = end;
            advance();
        }

        private void advance() {
            if (idx == features.size()) {
                next = null;
            } else {
                while (idx < features.size()) {
                    next = features.get(idx++);
                    if (next.getStart() > end) {
                        next = null;
                        break;   // Done
                    } else if (next.getEnd() > start) {
                        break;
                    }
                }
            }
        }

        @Override
        public boolean hasNext() {
            return next != null;
        }

        @Override
        public BasicFeature next() {
            BasicFeature retValue = next;
            advance();
            return retValue;
        }
    }
}
