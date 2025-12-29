package org.igv.feature.bionano;

import org.igv.logging.*;
import org.igv.feature.AbstractFeature;
import org.igv.track.WindowFunction;


/**
 * Created by jrobinson on 2/25/16.
 */
public class SMAPPairedFeature extends AbstractFeature {

    private static Logger log = LogManager.getLogger(SMAPPairedFeature.class);

    SMAPFeature feature1;
    SMAPFeature feature2;

    String type;

    public SMAPPairedFeature(SMAPFeature feature1, SMAPFeature feature2) {


        if (!feature1.getChr().equals(feature2.getChr())) {
            // TODO - throw error?
            log.error("Inter-chromosomal linked features not supported");
            return;
        }

        setChr(feature1.getChr());
        if (feature1.getStart() < feature2.getStart()) {
            this.feature1 = feature1;
            this.feature2  = feature2;
        } else {
            this.feature1 = feature2;
            this.feature2  = feature1;
        }

        // Feaures shouldn't overlap, but just in case check both
        setStart(Math.min(feature1.getStart(), feature2.getStart()));
        setEnd(Math.max(feature1.getEnd(), feature2.getEnd()));

    }


    @Override
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {

        StringBuffer buf = new StringBuffer();
        buf.append(feature1.getValueString(position, mouseX, windowFunction));
        buf.append("<hr>");
        buf.append(feature2.getValueString(position, mouseX, windowFunction));
        return buf.toString();

    }

}

