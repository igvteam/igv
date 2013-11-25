/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
