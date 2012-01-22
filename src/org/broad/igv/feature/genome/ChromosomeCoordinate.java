package org.broad.igv.feature.genome;

/**
* @author Jim Robinson
* @date 1/21/12
*/
public class ChromosomeCoordinate {

    private String chr;
    private int coordinate;

    public ChromosomeCoordinate(String chr, int coordinate) {
        this.chr = chr;
        this.coordinate = coordinate;
    }

    public String getChr() {
        return chr;
    }

    public int getCoordinate() {
        return coordinate;
    }
}
