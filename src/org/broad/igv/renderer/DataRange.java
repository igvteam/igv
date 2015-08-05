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
package org.broad.igv.renderer;

import org.broad.igv.track.Track;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import java.util.Collection;

/**
 * Encapsulates parameter for an x-y plot axis.
 *
 * @author jrobinso
 */
@XmlAccessorType(XmlAccessType.NONE)
public class DataRange{


    public enum Type {
        LOG, LINEAR
    }

    /**
     * The scale type,  linear by default
     */
    @XmlAttribute private Type type = Type.LINEAR;

    /**
     * Minimum data value displayed.  Zero by default.
     */
    @XmlAttribute private float minimum = 0;

    /**
     * Where to draw the plot baseline.  Zero by default
     */
    @XmlAttribute private float baseline = 0;

    /**
     * Maximum data value displayed. This value is required, no default
     */
    @XmlAttribute private float maximum;

    /**
     * If true the Y axis is "flipped" (most negative value at top)
     */
    @XmlAttribute private boolean flipAxis = false;

    @XmlAttribute private boolean drawBaseline = true;

    //Here for JAXB Compatibility
    private DataRange(){}

    public DataRange(float minimum, float maximum) {
        this(minimum, minimum, maximum, true);
    }


    public DataRange(float minimum, float baseline, float maximum) {
        this(minimum, baseline, maximum, true);
    }

    public DataRange(float minimum, float baseline, float maximum, boolean drawBaseline) {
        this(minimum, baseline, maximum, drawBaseline, false);
    }

    public DataRange(float minimum, float baseline, float maximum, boolean drawBaseline, boolean isLog) {
        this.minimum = minimum;
        this.baseline = baseline;
        this.maximum = maximum;
        this.drawBaseline = drawBaseline;
        this.type = isLog ? Type.LOG : Type.LINEAR;
    }

    public void setType(Type type) {
        this.type = type;
    }

    public Type getType() {
        return type;
    }

    public boolean isLog() {
        return type == Type.LOG;
    }


    public float getMinimum() {
        return minimum;
    }


    public float getBaseline() {
        return baseline;
    }


    public float getMaximum() {
        return maximum;
    }


    public boolean isFlipAxis() {
        return flipAxis;
    }

    public boolean isDrawBaseline() {
        return drawBaseline;
    }

    public void setDrawBaseline(boolean drawBaseline) {
        this.drawBaseline = drawBaseline;
    }

    public static DataRange getFromTracks(Collection<? extends Track> tracks){
        float min = Float.MAX_VALUE;
        float max = Float.MIN_VALUE;
        float mid = 0;
        boolean drawBaseline = true;
        boolean isLog = true;
        for (Track t : tracks) {
            DataRange dr = t.getDataRange();
            min = Math.min(min, dr.getMinimum());
            max = Math.max(max, dr.getMaximum());
            mid += dr.getBaseline();
            drawBaseline &= dr.isDrawBaseline();
            isLog &= dr.isLog();
        }
        mid /= tracks.size();
        if (mid < min) {
            mid = min;
        } else if (mid > max) {
            min = max;
        }

        return new DataRange(min, mid, max, drawBaseline, isLog);
    }

}
