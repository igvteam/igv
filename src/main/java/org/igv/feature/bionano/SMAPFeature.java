package org.igv.feature.bionano;

import org.igv.feature.AbstractFeature;
import org.igv.feature.Strand;
import org.igv.track.WindowFunction;

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
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {

        StringBuffer buf = new StringBuffer();
        buf.append("<b>Type:&nbsp;" + type + "</b>");
        for (int i = 0; i < headers.length; i++) {
            buf.append("<br>" + headers[i] + ":&nbsp;" + tokens[i]);
        }

        if(partialFeatures != null) {
            for(SMAPFeature pf : partialFeatures) {
                buf.append("<hr>");
                buf.append(pf.getValueString(position, mouseX, windowFunction));
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
