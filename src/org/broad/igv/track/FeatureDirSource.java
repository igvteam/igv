/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.feature.AbstractFeatureParser;
import org.broad.igv.feature.IGVFeature;
import org.broad.tribble.Feature;
import org.broad.igv.feature.FeatureParser;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
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

    public FeatureDirSource(ResourceLocator locator) throws IOException {
        featureCache = new LRUCache(this, 3);
        rootLocator = locator;
        setRootDir(locator.getPath());

        fileMap = new Properties();
        InputStream propStream = ParsingUtils.openInputStream(locator);
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
                AsciiLineReader reader = null;
                String path = rootDir + "/" + filename;
                try {
                    log.info("Loading " + path);
                    // Load features here
                    ResourceLocator loc = new ResourceLocator(rootLocator.getServerURL(), path);

                    FeatureParser fp = AbstractFeatureParser.getInstanceFor(loc);
                    reader = ParsingUtils.openAsciiReader(loc);
                    features = fp.loadFeatures(reader);
                    featureCache.put(chr, features);
                } catch (IOException ex) {
                    MessageUtils.showMessage("Error loading file: " + path + " (" + ex.toString() + ")");
                    log.info("Error loading feature file: " + filename, ex);
                } finally {
                    if (reader != null) {
                        reader.close();
                    }
                }

            }
        }
        return featureCache.get(chr);
    }

    public List<LocusScore> getCoverageScores(String chr, int i, int i1, int zoom) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getFeatureWindowSize() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setFeatureWindowSize(int size) {
        // ignored
    }


    public Iterator<Feature> getFeatures(String chr, int start, int end) {
        List<Feature> features = getFeatures(chr);
        return features == null ? null : features.iterator();
    }

    private void setRootDir(String path) {

        if (IGVHttpUtils.isURL(path)) {
            int idx = path.lastIndexOf('/');
            rootDir = path.substring(0, idx);
        } else {
            rootDir = (new File(path)).getParent();
        }

    }
}
