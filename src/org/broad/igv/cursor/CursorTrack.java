package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;

import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 *         Date: 1/14/14
 *         Time: 12:57 PM
 */
public class CursorTrack {


    Map<String, List<BasicFeature>> featureMap;
    Map<String, Integer> longestFeature;
    private Color color = new Color(0, 0, 150);
    private String name;
    Class featureType;

    public CursorTrack(Map<String, List<BasicFeature>> featureMap, Class featureType) {
        this.featureMap = featureMap;
        this.featureType = featureType;
        this.longestFeature = new HashMap();
        for(Map.Entry<String, List<BasicFeature>> entry : featureMap.entrySet()) {
            String chr = entry.getKey();
            List<BasicFeature> features = entry.getValue();
            int longest = 0;
            for(BasicFeature f : features) {
                final int length = f.getLength();
                if(length > longest) longest = length;
            }
            longestFeature.put(chr, longest);
        }
    }

    public List<BasicFeature> getFeatures(String chr) {
        return featureMap.get(chr);
    }

    public int getLongestFeatureLength(String chr) {
        Integer lf = longestFeature.get(chr);
        return lf == null ? -1 : lf;
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
