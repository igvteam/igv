package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * Wraps a collection of FeatureSource objects, one per chromosome.  The class is needed for 1KG variant files,
 * which are divided by chromosome.  It could be generalized to non-tribble sources, but no attempt has been
 * made to do so.
 *
 * @author jrobinso
 *         Date: 10/4/12
 *         Time: 9:36 PM
 */
public class TribbleListFeatureSource implements FeatureSource {

    private static Logger log = Logger.getLogger(TribbleListFeatureSource.class);

     Map<String, TribbleFeatureSource> featureSourceMap;
    int windowSize = 1000;
    Object header;
    Genome genome;

    public TribbleListFeatureSource(String path, Genome genome) throws IOException, TribbleIndexNotFoundException {

        this.genome = genome;
        init(path, genome);

    }

    private void init(String path, Genome genome) throws IOException, TribbleIndexNotFoundException {

        featureSourceMap = Collections.synchronizedMap(new HashMap());
        BufferedReader reader = null;

        try {
            reader = ParsingUtils.openBufferedReader(path);
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                String f = nextLine.trim();
                if (!f.startsWith("#")) {
                    String[] tokens = Globals.whitespacePattern.split(nextLine);
                    final String chr = tokens[0];
                    final String srcPath = tokens[1];
                    featureSourceMap.put(chr, TribbleFeatureSource.getFeatureSource(new ResourceLocator(srcPath), genome));

                }
            }
        } finally {
            if (reader != null) reader.close();
        }
    }



    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        FeatureSource src = featureSourceMap.get(chr);
        if (src != null) {
            return src.getFeatures(chr, start, end);
        } else {
            return null;
        }
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {

        FeatureSource src = featureSourceMap.get(chr);
        if (src != null) {
            return src.getCoverageScores(chr, start, end, zoom);
        } else {
            return null;
        }

    }

    @Override
    public int getFeatureWindowSize() {
        return windowSize;
    }

    @Override
    public void setFeatureWindowSize(int size) {
        this.windowSize = size;
    }

    public Object getHeader() {
        if (header == null) {
            // Arbitrarily get the first source
            if (featureSourceMap != null && featureSourceMap.size() > 0) {
                TribbleFeatureSource src = featureSourceMap.values().iterator().next();
                header = src.getHeader();
            }

        }
        return header;
    }
}
