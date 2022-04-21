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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class DotAlignedAlignment implements Alignment {

    String readName;
    private String chromosome;
    private int start;
    private int end;
    boolean negativeStrand;


    public DotAlignedAlignment(String chromosome, int start, int end, boolean isNegative, String name) {
        this.negativeStrand = isNegative;
        this.readName = name;
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
    }

    public DotAlignedAlignment(String chromosome, int start, int end, boolean negativeStrand) {
        this.readName = chromosome + ":" + (start + 1) + "-" + (end + 1) + "(" +
                (negativeStrand ? "-" : "+") + ")";
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.negativeStrand = negativeStrand;
    }

    public String getReadName() {
        return readName;
    }

    public boolean isSmallInsert() {
        return false;
    }

    public Color getYcColor() {
        return null;
    }


    public String getChromosome() {
        return chromosome;
    }

    public String getChr() {
        return chromosome;
    }

    @Override
    public String getContig() {
        return chromosome;
    }

    @Override
    public int getAlignmentStart() {
        return getStart();
    }

    public float getScore() {
        return 1.0f;
    }

    public LocusScore copy() {
        return this;
    }

    public String getClipboardString(double location, int mouseX) {
        return getValueString(location, mouseX, (WindowFunction) null);
    }

    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        return readName + "<br>Read length = " + (getEnd() - getStart());
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
        return end;
    }

    public int getAlignmentEnd() {
        return end;
    }

    /**
     * @param end the end to set
     */
    public void setEnd(int end) {
        this.end = end;
    }

    public Strand getReadStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }

}
