package org.broad.igv.ui.color;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by jrobinson on 7/5/16.
 */
public class GreyscaleColorTable implements ColorTable {

    Map<String, Color> colorMap;
    int min;
    int max;

    public GreyscaleColorTable() {
        this(50, 200);
    }

    public GreyscaleColorTable(int min, int max) {
        colorMap = new HashMap<>();
        this.min = min;
        this.max = max;
    }

    @Override
    public Color get(String key) {
        Color c = colorMap.get(key);
        if (c == null) {
            int r = (int) (min + Math.random() * (max - min));
            c = new Color(r, r, r);
            colorMap.put(key, c);

        }
        return c;
    }
}
