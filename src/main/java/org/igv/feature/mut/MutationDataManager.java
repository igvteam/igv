package org.igv.feature.mut;

import htsjdk.tribble.Feature;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.ParsingUtils;
import org.igv.Globals;
import org.igv.feature.Range;
import org.igv.feature.genome.Genome;
import org.igv.feature.tribble.MUTCodec;
import org.igv.feature.tribble.TribbleIndexNotFoundException;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.TribbleFeatureSource;
import org.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.*;

public class MutationDataManager {

    private static Logger log = LogManager.getLogger(MutationDataManager.class);

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

    private Index loadIndex(String path) {
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
