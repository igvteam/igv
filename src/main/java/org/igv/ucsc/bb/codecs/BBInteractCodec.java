package org.igv.ucsc.bb.codecs;

import org.igv.Globals;
import org.igv.bedpe.BedPEFeature;
import org.igv.bedpe.InteractFeature;
import org.igv.feature.BasicFeature;
import org.igv.feature.Exon;
import org.igv.feature.tribble.IGVBEDCodec;
import org.igv.ucsc.bb.BBUtils;
import org.igv.ucsc.bb.BedData;
import software.amazon.awssdk.services.s3.endpoints.internal.Value;

import java.util.List;


public class BBInteractCodec implements BBCodec {


    private final BBUtils.ASTable astable;

    private final int standardFieldCount;

    public BBInteractCodec(int standardFieldCount, BBUtils.ASTable autosql) {
        this.astable = autosql;
        this.standardFieldCount = standardFieldCount;
    }

    public InteractFeature decode(BedData bedData) {

        String[] tokens = bedData.toString().split("\t");

        return InteractFeature.fromTokens(tokens, null);

    }

}
