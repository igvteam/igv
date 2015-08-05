/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


/*
* GisticScore.java
*
* Created on June 20, 2007, 3:00 PM
*
* To change this template, choose Tools | Template Manager
* and open the template in the editor.
*/

package org.broad.igv.feature;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.track.WindowFunction;

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

    public String getValueString(double position, WindowFunction windowFunction) {
        return "";
    }


}
