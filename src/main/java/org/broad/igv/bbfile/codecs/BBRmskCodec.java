package org.broad.igv.bbfile.codecs;

import org.broad.igv.bbfile.BBUtils;
import org.broad.igv.bbfile.BedData;
import org.broad.igv.feature.BasicFeature;

import java.awt.*;
import java.util.Map;

public class BBRmskCodec implements BBCodec {

    BBBedCodec bedCodec;

    public BBRmskCodec(int standardFieldCount, BBUtils.ASTable autosql) {
        bedCodec = new BBBedCodec(standardFieldCount, autosql);
    }

    @Override
    public BasicFeature decode(BedData data) {
        BasicFeature feature = bedCodec.decode(data);

        String name = feature.getName();
        String repClass = null;
        String[] parts = name.split("#");
        if (parts.length > 1) {
            name = parts[0];
            repClass = parts[1];
        }

        feature.setName(name);
        if (repClass != null) {
            feature.setAttribute("Repeat Class", repClass);
            int idx1 = repClass.indexOf("/");
            String c = idx1 > 0 ? repClass.substring(0, idx1) : repClass;
            if (c.endsWith("?")) c = c.substring(c.length() - 1);
            Color color = classColors.containsKey(c) ? classColors.get(c) : defaultColor;
            feature.setColor(color);
        }

        return feature;
    }


    static Map<String, Color> classColors = Map.of(
            "SINE", new Color(31, 119, 180),
            "LINE", new Color(255, 127, 14),
            "LTR", new Color(44, 160, 44),
            "DNA", new Color(214, 39, 40),
            "Simple", new Color(148, 103, 189),
            "Low_complexity", new Color(140, 86, 75),
            "Satellite", new Color(227, 119, 194),
            "RNA", new Color(127, 127, 127),
            "Other", new Color(188, 189, 34),
            "Unknown", new Color(23, 190, 207)
    );

    static Color defaultColor = new Color(0, 180, 180);
}
