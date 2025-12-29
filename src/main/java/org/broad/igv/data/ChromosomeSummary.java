/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

/**
 * @author jrobinso
 */
public class ChromosomeSummary {

    /**
     * The chromosome name
     */
    private String name;
    /**
     * The approximate number of data points on this chromsome.  Used to estimate data
     * loading completion %.
     */
    private int nDataPoints;
    /**
     * The approximate number of characters to the starting position for this
     * chromosomes data in the copy number file.
     */
    private long startPosition;

    /**
     * Creates a new instance of ChromsomeSummary
     */
    public ChromosomeSummary(String name, long startPosition) {

        // Convert name to ucsc convention
        this.name = name; //.startsWith("chr") ? name : "chr" + name;
        this.startPosition = startPosition;
    }

    public String getName() {
        return name;
    }

    public long getStartPosition() {
        return startPosition;
    }

    public int getNDataPts() {
        return nDataPoints;
    }

    public void setNDataPoints(int nSnps) {
        this.nDataPoints = nSnps;
    }


}
