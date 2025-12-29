package org.igv.bedpe;

import org.igv.bedpe.InteractFeature;
import org.igv.feature.*;
import org.igv.feature.genome.Genome;
import org.igv.feature.tribble.UCSCCodec;

/**
 * Decode an UCSC interact file
 *
 * Reference: https://genome.ucsc.edu/goldenpath/help/interact.html
 *
 */
public class InteractCodec extends UCSCCodec<InteractFeature> {

    private Genome genome;

    public InteractCodec(Genome genome, FeatureType featureType) {
        super(InteractFeature.class, featureType);
        this.genome = genome;
    }


    //@Override
    public InteractFeature decode(String[] tokens) {
        return InteractFeature.fromTokens(tokens, genome);
    }

    public boolean canDecode(String path) {
        return true;
    }
}


