package org.broad.igv.bbfile.codecs;

import org.broad.igv.Globals;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;

import java.util.List;

public class BBGenePredCodec extends BBBedCodec {

    public BBGenePredCodec(int standardFieldCount) {
        super(standardFieldCount);
    }

    @Override
    public BasicFeature decode(BedFeature feat) {

        BasicFeature feature = super.decode(feat);

        String[] restOfFields = feat.getRestOfFields();
        String exonFrames = restOfFields[15 - 3];
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

        return feature;
    }
}
