
// should perhaps be in igv.data or some other location
package org.broad.igv.feature.basepair;

import org.broad.igv.feature.FeatureUtils;

import java.awt.*; // does this have Color definition?
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

public class BasePairData{

    Map<String, List<BasePairFeature>> featureMap;

    public BasePairData() {
        this.featureMap = new HashMap<String, List<BasePairFeature>>();
    }


    public void addFeature(BasePairFeature feature) {

        String chr = feature.getChr();
        List<BasePairFeature> featureList = featureMap.get(chr);
        if(featureList == null) {
            featureList = new ArrayList<BasePairFeature>();
            featureMap.put(chr, featureList);
        }
        featureList.add(feature);
    }

    public void finish() {
        for(List<BasePairFeature> featureList : featureMap.values()) {
            FeatureUtils.sortFeatureList(featureList);
        }
    }

    public List<BasePairFeature> getFeatures(String chr) {
        return featureMap.get(chr);
    }
}