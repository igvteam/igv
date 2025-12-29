package org.igv.feature.dsi;

import htsjdk.tribble.Feature;
import org.igv.track.FeatureSource;
import org.igv.track.FeatureTrack;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;

import java.util.List;

/**
 * Created by jrobinson on 7/19/16.
 */
public class DSITrack extends FeatureTrack {


    public DSITrack() {
    }

    public DSITrack(ResourceLocator locator, FeatureSource src) {
        super(locator, src);
        setRenderer(new DSIRenderer());
    }


    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        List<Feature> allFeatures = getAllFeaturesContaining(position, mouseY, frame);
        if (allFeatures == null) {
            return null;
        }

        StringBuffer buf = new StringBuffer();
        boolean firstFeature = true;
        int maxNumber = 10;
        int n = 1;
        for (Feature feature : allFeatures) {

            if (feature != null && feature instanceof DSIFeature) {

                if (!firstFeature) {
                    buf.append("<hr><br>");
                }

                DSIFeature igvFeature = (DSIFeature) feature;
                String vs = igvFeature.getValueString(position, null);
                buf.append(vs);


                firstFeature = false;

                if (n > maxNumber) {
                    buf.append("...");
                    break;
                }
            }
            n++;
        }

        return buf.toString();
    }
}
