package org.broad.igv.feature.gff;

import org.broad.igv.logging.*;
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

    private static Logger log = LogManager.getLogger(GFFFeatureSource.class);
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
    public Object getHeader() throws IOException {
        return this.wrappedSource.getHeader();
    }

    @Override
    public void close() {
        wrappedSource.close();
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

}
