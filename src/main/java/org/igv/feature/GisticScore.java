/*
* GisticScore.java
*
* Created on June 20, 2007, 3:00 PM
*
* To change this template, choose Tools | Template Manager
* and open the template in the editor.
*/

package org.igv.feature;


import org.igv.track.WindowFunction;

/**
 * @author jrobinso
 */
public class GisticScore implements LocusScore {

    public enum Type {
        AMP, DEL
    }

    ;

    // The chromosome/start/end triplet could be encapulated as a "location"
    private String chromosome;
    int start;
    int end;
    private float qValue;
    private float gScore;
    private Type gisticType;

    /**
     * Creates a new instance of GisticScore
     */
    public GisticScore() {
    }

    /**
     * Constructs ...
     *
     * @param chr
     * @param start
     * @param end
     * @param qValue
     * @param gScore
     * @param type
     */
    public GisticScore(String chr, int start, int end, float qValue,
                       float gScore, Type type) {
        this.chromosome = chr;
        this.start = start;
        this.end = end;
        this.qValue = qValue;
        this.gScore = gScore;
        this.gisticType = type;
    }

    public GisticScore(GisticScore gisticScore) {
        this.chromosome = gisticScore.chromosome;
        this.start = gisticScore.start;
        this.end = gisticScore.end;
        this.qValue = gisticScore.qValue;
        this.gScore = gisticScore.gScore;
        this.gisticType = gisticScore.gisticType;
    }

    public GisticScore copy() {
        return new GisticScore(this);
    }

    public String getChromosome() {
        return chromosome;
    }

    public float getScore() {
        return gScore;
    }

    public double getQValue() {
        return qValue;
    }

    public double getGScore() {
        return gScore;
    }

    public void setQValue(float qValue) {
        this.qValue = qValue;
    }

    public void setGScore(float gScore) {
        this.gScore = gScore;
    }

    public Type getType() {
        return gisticType;
    }

    public void setType(Type type) {
        this.gisticType = type;
    }

    public String getChr() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getContig() {
        return null;
    }

    /**
     * Method description
     *
     * @return
     */
    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }


    public void setEnd(int end) {
        this.end = end;
    }

    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        return "";
    }


}
