package org.broad.igv.ucsc.bb.codecs;

import org.broad.igv.Globals;
import org.broad.igv.bedpe.BedPEFeature;
import org.broad.igv.bedpe.InteractFeature;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.ucsc.bb.BBUtils;
import org.broad.igv.ucsc.bb.BedData;
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
