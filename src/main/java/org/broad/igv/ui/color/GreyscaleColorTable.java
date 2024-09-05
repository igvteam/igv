package org.broad.igv.ui.color;

import java.awt.*;

/**
 * Created by jrobinson on 7/5/16.
 */
public class GreyscaleColorTable extends ColorTable {

    final int min;
    final int max;

    public GreyscaleColorTable() {
        this(50, 200);
    }

    public GreyscaleColorTable(int min, int max) {
        this.min = min;
        this.max = max;
    }

    @Override
    protected Color computeColor(String key) {
        int r = (int) (min + Math.random() * (max - min));
        return new Color(r, r, r);
    }
}
