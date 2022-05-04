package org.broad.igv.bbfile.codecs;

import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.tribble.IGVBEDCodec;

public class BBBedCodec implements BBCodec {

    int standardFieldCount;
    IGVBEDCodec igvBedCodec;

    public BBBedCodec(int standardFieldCount) {
        this.standardFieldCount = standardFieldCount;
        this.igvBedCodec = new IGVBEDCodec();   // Backing "tribble" codec
    }

    public BasicFeature decode(BedFeature feat) {

        String[] restOfFields = feat.getRestOfFields();
        String[] tokens = new String[restOfFields.length + 3];
        tokens[0] = feat.getChromosome();
        tokens[1] = String.valueOf(feat.getStartBase());
        tokens[2] = String.valueOf(feat.getEndBase());

        System.arraycopy(restOfFields, 0, tokens, 3, this.standardFieldCount - 3);
        BasicFeature feature = igvBedCodec.decode(tokens);
        return feature;

    }
}
