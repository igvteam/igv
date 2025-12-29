package org.igv.feature;

import org.igv.track.FeatureSource;
import htsjdk.tribble.Feature;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * A FeatureSource wrapper which provides caching.
 * The cache is only cleared when the source is closed.
 *
 * @author jacob
 * @date 2012/05/15
 */
public class CachingFeatureSource extends AbstractCacher implements FeatureSource {

    private static final int maxBinCount = 1000;
    private static final int defaultBinSize = 16000; // <= 16 kb

    private FeatureSource source;


    /**
     * Wraps the provided {@code source} with a caching version,
     * using default parameters.
     * @param source
     * @api
     */
    public CachingFeatureSource(FeatureSource source) {
        this(source, maxBinCount, defaultBinSize);
    }


    public CachingFeatureSource(FeatureSource source, int binCount, int binSize) {
        super(binCount, binSize);
        this.source = source;
    }

    @Override
    protected Iterator<Feature> queryRaw(String chr, int start, int end) throws IOException {
        return source.getFeatures(chr, start, end);
    }

    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        return super.queryCached(chr, start, end);
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return source.getCoverageScores(chr, start, end, zoom);
    }

    @Override
    public int getFeatureWindowSize() {
        return source.getFeatureWindowSize();
    }

    /**
     * Return the source which backs this CachingFeatureSource
     * @return
     */
    public FeatureSource getSource(){
        return this.source;
    }
}
