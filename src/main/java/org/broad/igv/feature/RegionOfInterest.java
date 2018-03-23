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


package org.broad.igv.feature;



import org.broad.igv.ui.WaitCursorManager;

import java.awt.*;

/**
 * @author eflakes
 */
public class RegionOfInterest{

    private String chr;
    private String description;
    private int start;    // In Chromosome coordinates
    private int end;      // In Chromosome coordinates
    private static Color backgroundColor = Color.RED;
    private static Color foregroundColor = Color.BLACK;
    boolean selected = false;

    private WaitCursorManager.CursorToken token;

    /**
     * A bounded region on a chromosome.
     *
     * @param chromosomeName
     * @param start          The region starting position on the chromosome.
     * @param end            The region starting position on the chromosome.
     * @param description
     */
    public RegionOfInterest(String chromosomeName, int start, int end, String description) {

        this.chr = chromosomeName;
        this.description = description;
        this.start = start;
        this.end = end;
    }

    public String getTooltip() {
        return description == null ? chr + ":" + getDisplayStart() + "-" + getDisplayEnd() : description;
    }

    public String getChr() {
        return chr;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getDescription() {
        return description;
    }


    public void setEnd(int end) {
        this.end = end;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    /**
     * locations displayed to the user are 1-based.  start and end are 0-based.
     * @return
     */
    public int getDisplayEnd() {
        return getEnd();
    }

    public int getStart() {
        return start;
    }

    public int getCenter() {
        return (start + end) / 2;
    }

    public int getLength() {
        return end - start;
    }

    /**
     * locations displayed to the user are 1-based.  start and end are 0-based.
     * @return
     */
    public int getDisplayStart() {
        return getStart() + 1;
    }

    public static Color getBackgroundColor() {
        return backgroundColor;
    }

    public static Color getForegroundColor() {
        return foregroundColor;
    }


    public String getLocusString() {
        return getChr() + ":" + getDisplayStart() + "-" + getDisplayEnd();
    }
}
