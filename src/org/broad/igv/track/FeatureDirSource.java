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
import org.broad.tribble.Feature;

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
        featureCache = new LRUCache(this, 3);
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
                    ResourceLocator loc = new ResourceLocator(rootLocator.getServerURL(), path);

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
