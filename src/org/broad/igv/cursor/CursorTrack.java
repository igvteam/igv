package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;

import java.awt.*;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 *         Date: 1/14/14
 *         Time: 12:57 PM
 */
public class CursorTrack {


    Map<String, List<BasicFeature>> featureMap;
    private Color color = new Color(0, 0, 150);
    private String name;
    Class featureType;

    public CursorTrack(Map<String, List<BasicFeature>> featureMap, Class featureType) {
        this.featureMap = featureMap;
        this.featureType = featureType;
    }

    public List<BasicFeature> getFeatures(String chr) {
        return featureMap.get(chr);
    }

    public Map<String, List<BasicFeature>> getFeatureMap() {
        return featureMap;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }


}
