package org.igv.feature;

import org.igv.feature.genome.Genome;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;

/**
 * FeatureParser which reads features using a codec.
 * Intended as a bridge between FeatureParser interface
 * and FeatureSource interface. Hope to remove the FeatureParser
 * interface eventually, making this no longer necessary.
 * <p/>
 * <p/>
 * <p/>
 * User: jacob
 * Date: 2012/01/26
 */
public class FeatureCodecParser extends AbstractFeatureParser {
    private AsciiFeatureCodec codec;

    public FeatureCodecParser(AsciiFeatureCodec codec, Genome genome) {
        this.codec = codec;
    }

    @Override
    protected Feature parseLine(String nextLine) {
        return codec.decode(nextLine);
    }
}
