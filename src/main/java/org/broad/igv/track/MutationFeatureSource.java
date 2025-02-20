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

package org.broad.igv.track;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.ParsingUtils;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Mutation;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.MUTCodec;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 *         Date: 4/9/13
 *         Time: 8:45 AM
 */
public class MutationFeatureSource implements FeatureSource<Mutation> {

    private static Logger log = LogManager.getLogger(MutationFeatureSource.class);

    String sample;
    MutationDataManager dataManager;


    public MutationFeatureSource(String sample, MutationDataManager dataManager) {
        this.sample = sample;
        this.dataManager = dataManager;
    }

    @Override
    public Iterator<Mutation> getFeatures(String chr, int start, int end) throws IOException {
        return dataManager.getFeatures(sample, chr, start, end);
    }

    static public class MutationDataManager {

        Range currentRange;
        Map<String, List<Mutation>> featureMap = Collections.synchronizedMap(new HashMap());
        TribbleFeatureSource tribbleFeatureSource;
        ResourceLocator locator;

        public MutationDataManager(ResourceLocator locator, Genome genome) throws IOException, TribbleIndexNotFoundException {
            this.tribbleFeatureSource = TribbleFeatureSource.getFeatureSource(locator, genome);
            this.locator = locator;
        }

        public boolean isIndexed() {
            return this.tribbleFeatureSource.isIndexed();
        }

        synchronized Iterator<Mutation> getFeatures(String trackKey, String chr, int start, int end) throws IOException {
            if (currentRange == null || !currentRange.contains(chr, start, end)) {
                Iterator<Feature> features = tribbleFeatureSource.getFeatures(chr, start, end);

                while (features.hasNext()) {
                    Mutation feat = (Mutation) features.next();
                    String thisKey = feat.getSampleId();
                    List<Mutation> keyFeatures = featureMap.get(thisKey);
                    if (keyFeatures == null) {
                        keyFeatures = new ArrayList<Mutation>();
                        featureMap.put(thisKey, keyFeatures);
                    }
                    keyFeatures.add(feat);
                    currentRange = new Range(chr, start, end);
                }

            }
            List<Mutation> featureList = featureMap.get(trackKey);
            return featureList == null ? Collections.EMPTY_LIST.iterator() : featureList.iterator();

        }


        public String[] getSamples() {

            if (this.tribbleFeatureSource.isIndexed()) {
                // Load the index for meta information.  Its already loaded, but not public in the htsjdk class.
                Index idx = loadIndex(this.locator.getPath());
                if (idx != null) {
                    Map<String, String> map = idx.getProperties();
                    if (map != null && map.containsKey("samples")) {
                        return Globals.commaPattern.split(map.get("samples"));
                    }
                }
            }

            // Try to fetch features from codec.  This is to support a deprecated option to
            // specify sample names in the .mut or .maf file header.

            MUTCodec codec = new MUTCodec(locator.getPath(), null);
            return codec.getSamples();

        }

        private Index loadIndex(String path)  {
            try {
                String indexFile = Tribble.indexFile(path);
                Index index = null;
                if (ParsingUtils.resourceExists(indexFile)) {
                    index = IndexFactory.loadIndex(indexFile);
                } else {
                    // See if the index itself is gzipped
                    indexFile = ParsingUtils.appendToPath(indexFile, ".gz");
                    if (ParsingUtils.resourceExists(indexFile)) {
                        index = IndexFactory.loadIndex(indexFile);
                    }
                }
                return index;
            } catch (IOException e) {
                log.error("Error loading index file", e);
                return null;
            }
        }
    }
}
