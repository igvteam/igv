/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
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

package org.broad.igv.feature.bionano;


import org.broad.igv.feature.AbstractFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class SMAPFeature extends AbstractFeature {

    private int linkId;
    private String[] headers;
    private String[] tokens;
    double confidence;
    List<SMAPFeature> partialFeatures;

    public SMAPFeature(String chr, int start, int end, double confidence, String type, String[] headers, String[] tokens) {
        super(chr, start, end, Strand.NONE);

        // TODO -- check tokens length
        this.tokens = tokens;
        this.confidence = confidence;
        this.type = type;
        this.headers = headers;

    }

    public SMAPFeature(String chr, int start, int end, double conf, String t, String[] headers, String[] tokens, int linkId) {
        this(chr, start, end, conf, t, headers, tokens);
        this.linkId = linkId;
    }

    @Override
    public Color getColor() {
        if (colors.containsKey(type)) {
            return colors.get(type);
        } else {
            return super.getColor();
        }
    }

    public void addPartialFeature(SMAPFeature smapFeature) {
        if (partialFeatures == null) {
            partialFeatures = new ArrayList<SMAPFeature>();
        }
        partialFeatures.add(smapFeature);
    }

    @Override
    public String getURL() {
        return null;
    }

    @Override
    public String getValueString(double position, WindowFunction windowFunction) {

        StringBuffer buf = new StringBuffer();
        buf.append("<b>Type:&nbsp;" + type + "</b>");
        for (int i = 0; i < headers.length; i++) {
            buf.append("<br>" + headers[i] + ":&nbsp;" + tokens[i]);
        }

        if(partialFeatures != null) {
            for(SMAPFeature pf : partialFeatures) {
                buf.append("<hr>");
                buf.append(pf.getValueString(position, windowFunction));
            }
        }

        return buf.toString();
    }

    public int getLinkId() {
        return linkId;
    }

    static Map<String, Color> colors = new HashMap<String, Color>();

    static {

        colors.put("insertion", new Color(0, 128, 0));
        colors.put("deletion", new Color(255, 0, 0));

    }
}
