package org.igv.feature;

import org.igv.track.WindowFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 * Date: 7/18/13
 * Time: 9:50 PM
 */
public class GFFFeature extends BasicFeature {

    List<String> componentAttributes = new ArrayList<>();

    public GFFFeature(BasicFeature feature) {
        super(feature);
    }

    public GFFFeature(String chr, int start, int end, Strand strand) {
        super(chr, start, end, strand);
    }

    @Override
    public void addExon(Exon region) {
        super.addExon(region);
    }

    @Override
    public void addUTRorCDS(BasicFeature bf) {
        super.addUTRorCDS(bf);
    }

    @Override
    public String getValueString(double position, int mouseX, WindowFunction ignored) {

        StringBuffer valueString = new StringBuffer();

        valueString.append("<b>type:</b>&nbsp;" + this.type);

        if (attributes != null) {
            valueString.append(getAttributeString());
        }

        if (componentAttributes.size() > 0) {
            for (String s : componentAttributes) {
                valueString.append("---------------------------");
                valueString.append(s);
            }
        }
        return valueString.toString();
    }

    public void mergeAttributes(BasicFeature mrnaPart) {
        StringBuffer buf = new StringBuffer();
        buf.append("<br><b>type:</b>&nbsp;" + mrnaPart.getType());
        FormatUtils.printAttributes(mrnaPart.getAttributes(), buf, 100);
        componentAttributes.add(buf.toString());
    }
}
