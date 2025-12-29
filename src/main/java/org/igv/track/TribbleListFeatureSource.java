package org.igv.track;

import org.igv.logging.*;
import org.igv.Globals;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.feature.tribble.TribbleIndexNotFoundException;
import org.igv.util.ParsingUtils;
import org.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.FileReader;
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

    private static Logger log = LogManager.getLogger(TribbleListFeatureSource.class);

    Map<String, String> pathMap;
    Map<String, TribbleFeatureSource> featureSourceMap;
    int windowSize = 1000;
    Object header;
    Genome genome;

    public TribbleListFeatureSource(String path, Genome genome) throws IOException {

        this.genome = genome;
        init(path, genome);

    }

    private void init(String path, Genome genome) throws IOException {

        featureSourceMap = Collections.synchronizedMap(new HashMap());
        pathMap = new HashMap<String, String>();
        BufferedReader reader = null;

        try {
            reader = ParsingUtils.openBufferedReader(path);
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                String f = nextLine.trim();
                if (!f.startsWith("#")) {
                    String[] tokens = Globals.whitespacePattern.split(nextLine);
                    if (tokens.length > 1)
                        pathMap.put(tokens[0], tokens[1]);
                }
            }
        } finally {
            if (reader != null) reader.close();
        }
    }

    private TribbleFeatureSource getSource(String chr) throws TribbleIndexNotFoundException {

        TribbleFeatureSource src = featureSourceMap.get(chr);
        if (src == null) {
            String path = pathMap.get(chr);
            if (path != null) {
                try {
                    src = TribbleFeatureSource.getFeatureSource(new ResourceLocator(path), genome);
                } catch (IOException e) {
                    log.error("Error loading tribble source: " + path);
                }
                featureSourceMap.put(chr, src);
            }
        }
        return src;

    }

    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        try {
            FeatureSource src = getSource(chr);
            if (src != null) {
                return src.getFeatures(chr, start, end);
            } else {
                return null;
            }
        } catch (TribbleIndexNotFoundException e) {
            log.error("Index not found", e);
            throw new IOException(e.getMessage());
        }

    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {

        try {
            FeatureSource src = getSource(chr);
            if (src != null) {
                return src.getCoverageScores(chr, start, end, zoom);
            } else {
                return null;
            }
        } catch (TribbleIndexNotFoundException e) {
            log.error("Index not found for", e);
            return null;
        }

    }

    @Override
    public int getFeatureWindowSize() {
        return windowSize;
    }

    public Object getHeader()  {
        if (header == null) {
            // Arbitrarily get the first source
            if (pathMap != null && pathMap.size() > 0) {
                String chr = pathMap.keySet().iterator().next();
                TribbleFeatureSource src = null;
                try {
                    src = getSource(chr);
                } catch (TribbleIndexNotFoundException e) {
                    log.error("Error reading header form Tribble source", e);
                    return null;
                }
                header = src.getHeader();
            }
        }
        return header;
    }
}
