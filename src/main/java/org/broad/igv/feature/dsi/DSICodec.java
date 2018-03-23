package org.broad.igv.feature.dsi;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.exception.CodecLineParsingException;
import htsjdk.tribble.readers.LineIterator;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.ParsingUtils;


/**
 * Created by jrobinson on 7/19/16.
 */
public class DSICodec extends AsciiFeatureCodec<DSIFeature> {

    private final Genome genome;
    private TrackProperties trackProperties;


    public DSICodec(Genome genome) {
        super(DSIFeature.class);
        this.genome = genome;
    }

    @Override
    public DSIFeature decode(String line) {

        String [] tokens = Globals.tabPattern.split(line);
        if(tokens.length > 10) {

            String chr = genome == null ? tokens[0].trim() : genome.getCanonicalChrName(tokens[0].trim());

            DSIFeature feature = new DSIFeature();
            feature.chr = chr;
            feature.position = Integer.parseInt(tokens[1]) - 1;   // IGV is zero-based internally
            feature.base = tokens[2].charAt(0);
            feature.total = Integer.parseInt(tokens[3]);
            feature.meth = Integer.parseInt(tokens[4]);
            feature.unmeth = Integer.parseInt(tokens[5]);
            feature.type = tokens[6];
            feature.f = "NA".equals(tokens[7]) ? Integer.MIN_VALUE : Integer.parseInt(tokens[7]);
            feature.p = "NA".equals(tokens[7]) ? Integer.MIN_VALUE : Integer.parseInt(tokens[8]);
            feature.m = "NA".equals(tokens[7]) ? Integer.MIN_VALUE : Integer.parseInt(tokens[9]);
            feature.u = "NA".equals(tokens[7]) ? Integer.MIN_VALUE : Integer.parseInt(tokens[10]);
            return feature;
        }
        else {
            return null;
        }
    }


    @Override
    public boolean canDecode(String path) {
        return path.toLowerCase().endsWith(".dsi") || path.toLowerCase().endsWith(".dsi.gz");

    }

    @Override
    public Object readActualHeader(LineIterator reader) {

        String line;
        try {
            while (reader.hasNext()) {
                line = reader.peek();
                if (line.startsWith("#")) {
                    reader.next();
                } else if (line.startsWith("#track") || line.startsWith("##track")) {
                    trackProperties = new TrackProperties();
                    ParsingUtils.parseTrackLine(line, trackProperties);
                } else {
                    break;
                }
            }

            return trackProperties;
        } catch (Exception e) {
            throw new CodecLineParsingException("Error parsing header: " + e.getMessage(), e);
        }
    }
}
