package org.igv.feature;

import htsjdk.tribble.Feature;

public interface PackedFeature extends Feature {

    void setPackedRow(int rowIndex);

    int getPackedRow();
}
