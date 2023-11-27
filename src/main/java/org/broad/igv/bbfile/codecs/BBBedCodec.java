package org.broad.igv.bbfile.codecs;

import org.broad.igv.Globals;
import org.broad.igv.ucsc.bb.BBUtils;
import org.broad.igv.bbfile.BedData;
import org.broad.igv.feature.*;
import org.broad.igv.feature.tribble.IGVBEDCodec;

import java.util.LinkedHashMap;
import java.util.List;


public class BBBedCodec implements BBCodec {


    private final BBUtils.ASTable astable;

    int standardFieldCount;
    IGVBEDCodec igvBedCodec;

    public BBBedCodec(int standardFieldCount, BBUtils.ASTable autosql) {
        this.astable = autosql;
        this.standardFieldCount = standardFieldCount;
        this.igvBedCodec = new IGVBEDCodec();   // Backing "tribble" codec
    }

    public BasicFeature decode(BedData bedData) {

        String[] restOfFields = Globals.tabPattern.split(bedData.getRestOfFields(), -1);
        String[] tokens = new String[this.standardFieldCount];
        tokens[0] = bedData.getChromosome();
        tokens[1] = String.valueOf(bedData.getStartBase());
        tokens[2] = String.valueOf(bedData.getEndBase());

        System.arraycopy(restOfFields, 0, tokens, 3, this.standardFieldCount - 3);
        BasicFeature feature = igvBedCodec.decode(tokens);

        // "Non standard" fields
        if (astable != null && astable.fields.size() > standardFieldCount) {
            LinkedHashMap<String, String> attributes = new LinkedHashMap<>();
            List<BBUtils.ASField> fields = astable.fields;
            for (int i = this.standardFieldCount; i < fields.size(); i++) {
                BBUtils.ASField field = fields.get(i);
                attributes.put(field.name, restOfFields[i - 3]);
            }
            feature.setAttributes(attributes);

            if (attributes.containsKey("exonFrames")) {
                computeExonFrames(feature, attributes.get("exonFrames"));
            }
        }

        return feature;
    }

    private void computeExonFrames(BasicFeature feature, String exonFrames) {
        String[] frameBuffer = Globals.commaPattern.split(exonFrames);
        List<Exon> exons = feature.getExons();
        for (int i = 0; i < frameBuffer.length; i++) {
            int exonFrame = Integer.parseInt(frameBuffer[i].trim());
            if (exonFrame == -1) {
                exons.get(i).setNonCoding(true);
            } else {
                exons.get(i).setReadingFrame(exonFrame);
            }
        }
    }


}
