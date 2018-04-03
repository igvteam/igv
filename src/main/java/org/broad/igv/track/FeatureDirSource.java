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

import org.apache.log4j.Logger;
import org.broad.igv.feature.AbstractFeatureParser;
import org.broad.igv.feature.FeatureParser;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.LRUCache;
import htsjdk.tribble.Feature;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

/**
 * User: jrobinso
 * Date: Jan 31, 2010
 */
public class FeatureDirSource implements FeatureSource {

    static Logger log = Logger.getLogger(FeatureDirSource.class);
    LRUCache<String, List<Feature>> featureCache;
    Properties fileMap;
    String rootDir;
    ResourceLocator rootLocator;
    Genome genome;

    public FeatureDirSource(ResourceLocator locator, Genome genome) throws IOException {
        this.genome = genome;
        featureCache = new LRUCache(3);
        rootLocator = locator;
        setRootDir(locator.getPath());

        fileMap = new Properties();
        InputStream propStream = ParsingUtils.openInputStreamGZ(locator);
        fileMap.load(propStream);
        propStream.close();


    }

    public Class getFeatureClass() {
        return IGVFeature.class;
    }


    public List<Feature> getFeatures(final String chr) {

        List<Feature> features = featureCache.get(chr);
        if (features == null) {
            final String filename = fileMap.getProperty(chr);
            if (filename != null) {
                BufferedReader reader = null;
                String path = rootDir + "/" + filename;
                try {
                    log.info("Loading " + path);
                    // Load features here
                    ResourceLocator loc = new ResourceLocator(path);

                    FeatureParser fp = AbstractFeatureParser.getInstanceFor(loc, genome);
                    reader = ParsingUtils.openBufferedReader(loc);
                    features = fp.loadFeatures(reader, genome);
                    featureCache.put(chr, features);
                } catch (IOException ex) {
                    MessageUtils.showMessage("Error loading file: " + path + " (" + ex.toString() + ")");
                    log.info("Error loading feature file: " + filename, ex);
                } finally {
                    if (reader != null) {
                        try {
                            reader.close();
                        } catch (IOException e) {

                        }
                    }
                }

            }
        }
        return featureCache.get(chr);
    }

    public List<LocusScore> getCoverageScores(String chr, int i, int i1, int zoom) {
        return null;
    }

    public int getFeatureWindowSize() {
        return 0;
    }

    public void setFeatureWindowSize(int size) {
        // ignored
    }


    public Iterator<Feature> getFeatures(String chr, int start, int end) {
        List<Feature> features = getFeatures(chr);
        return features == null ? Collections.<Feature>emptyList().iterator() : features.iterator();
    }

    private void setRootDir(String path) {

        if (HttpUtils.isRemoteURL(path)) {
            int idx = path.lastIndexOf('/');
            rootDir = path.substring(0, idx);
        } else {
            rootDir = (new File(path)).getParent();
        }

    }
}
