package org.igv.track;

import org.igv.logging.*;
import org.igv.feature.AbstractFeatureParser;
import org.igv.feature.FeatureParser;
import org.igv.feature.IGVFeature;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.ui.util.MessageUtils;
import org.igv.util.HttpUtils;
import org.igv.util.ParsingUtils;
import org.igv.util.ResourceLocator;
import org.igv.util.collections.LRUCache;
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

    static Logger log = LogManager.getLogger(FeatureDirSource.class);
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
                    //log.warn("Loading " + path);
                    // Load features here
                    ResourceLocator loc = new ResourceLocator(path);

                    FeatureParser fp = AbstractFeatureParser.getInstanceFor(loc, genome);
                    reader = ParsingUtils.openBufferedReader(loc);
                    features = fp.loadFeatures(reader, genome);
                    featureCache.put(chr, features);
                } catch (IOException ex) {
                    MessageUtils.showMessage("Error loading file: " + path + " (" + ex.toString() + ")");
                    log.error("Error loading feature file: " + filename, ex);
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
