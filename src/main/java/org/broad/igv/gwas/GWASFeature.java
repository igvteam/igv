package org.broad.igv.gwas;

import htsjdk.tribble.Feature;

public class GWASFeature  {


    String chr;
    int position;
    double value;
    String line;
    int pixelX;
    int pixelY;

    public GWASFeature(String chr, int position, double value, String line) {
        this.chr = chr;
        this.position = position;
        this.value = value;
        this.line = line;
    }


}
