/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.cursor;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.EncodePeakFeature;
import org.broad.igv.feature.SignalFeature;
import org.broad.igv.util.collections.DoubleArrayList;

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


    public enum SignalField {Score, Signal, PValue, QValue}


    Map<String, List<BasicFeature>> featureMap;
    Map<String, Integer> longestFeature;
    private Color color = new Color(0, 0, 150);
    private String name;
    Class featureType;
    boolean hasScore = false;
    private SignalField signalField = SignalField.Score;
    private Map<SignalField, Range> scales;

    public static class Range {
        double min;
        double max;

        private Range(double min, double max) {
            this.min = min;
            this.max = max;
        }

        public double getMax() {
            return max;
        }

        public double getMin() {
            return min;
        }
    }

    public CursorTrack(Map<String, List<BasicFeature>> featureMap, Class featureType) {
        this.featureMap = featureMap;
        this.featureType = featureType;
        this.longestFeature = new HashMap();

        // Gather statistics for setting scales.
        boolean isEncodePeak = featureType == EncodePeakFeature.class;
        DoubleArrayList signals = isEncodePeak ? new DoubleArrayList(50000) : null;
        DoubleArrayList qvalues = isEncodePeak ? new DoubleArrayList(50000) : null;
        DoubleArrayList pvalues = isEncodePeak ? new DoubleArrayList(50000) : null;

        for (Map.Entry<String, List<BasicFeature>> entry : featureMap.entrySet()) {
            String chr = entry.getKey();
            List<BasicFeature> features = entry.getValue();
            int longest = 0;
            for (BasicFeature f : features) {

                float s = f.getScore();
                hasScore = hasScore || !Float.isNaN(s);
                if (isEncodePeak) {
                    EncodePeakFeature pf = (EncodePeakFeature) f;
                    signals.add(pf.getSignal());
                    qvalues.add(pf.getQValue());
                    pvalues.add(pf.getPValue());
                }

                final int length = f.getLength();
                if (length > longest) longest = length;


            }
            longestFeature.put(chr, longest);
        }

        scales = new HashMap<SignalField, Range>();
        scales.put(SignalField.Score, new Range(0, 1000));

        double[] s = signals.toArray();
        double sMin = 0;
        double sMax = StatUtils.percentile(s, 90);
        scales.put(SignalField.Signal, new Range(sMin, sMax));

        double[] q = qvalues.toArray();
        scales.put(SignalField.QValue, new Range(0, StatUtils.percentile(q, 90)));

        double[] p = pvalues.toArray();
        scales.put(SignalField.PValue, new Range(0, StatUtils.percentile(p, 90)));

    }


    public Range getScale() {
        return scales.get(signalField);
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

    public SignalField getSignalField() {
        return signalField;
    }

    public void setSignalField(SignalField signalField) {
        System.out.println("Set signal field: " + signalField);
        this.signalField = signalField;
    }

    public SignalField[] getAllSignalFields() {
        if (!hasScore) return noSignalsArray;
        else if (featureType == EncodePeakFeature.class) return signalsArray;
        else return scoreOnlyArray;
    }


    public float getSignal(BasicFeature feature) {

        if (signalField == SignalField.Score) {
            return feature.getScore();
        } else {
            SignalFeature sf = (SignalFeature) feature;
            switch (signalField) {
                case Signal:
                    return sf.getSignal();
                case QValue:
                    return sf.getQValue();
                case PValue:
                    return sf.getPValue();
                default:
                    return feature.getScore();

            }
        }
    }

    // Some static arrays for effecienty
    private static SignalField[] noSignalsArray = new SignalField[0];
    private static SignalField[] scoreOnlyArray = new SignalField[]{SignalField.Score};
    private static SignalField[] signalsArray = new SignalField[]{
            SignalField.Score, SignalField.Signal, SignalField.PValue, SignalField.QValue};
}
