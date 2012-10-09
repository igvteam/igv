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

import org.broad.igv.session.Persistable;

import java.util.HashMap;
import java.util.Map;

/**
 * Encapsulates parameter for an x-y plot axis.
 *
 * @author jrobinso
 */
public class DataRange implements Persistable {

    public enum Type {
        LOG, LINEAR
    }

    /**
     * The scale type,  linear by default
     */
    private Type type = Type.LINEAR;

    /**
     * Minimum data value displayed.  Zero by default.
     */
    private float minimum = 0;

    /**
     * Where to draw the plot baseline.  Zero by default
     */
    private float baseline = 0;

    /**
     * Maximum data value displayed. This value is required, no default
     */
    private float maximum;

    /**
     * If true the Y axis is "flipped" (most negative value at top)
     */
    private boolean flipAxis = false;

    private boolean drawBaseline = true;


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


    public Map<String, String> getPersistentState() {
        Map<String, String> attributes = new HashMap();
        attributes.put("type", type.toString());
        attributes.put("minimum", String.valueOf(minimum));
        attributes.put("baseline", String.valueOf(baseline));
        attributes.put("maximum", String.valueOf(maximum));
        attributes.put("flipAxis", String.valueOf(flipAxis));
        attributes.put("drawBaseline", String.valueOf(drawBaseline));
        return attributes;
    }

    public void restorePersistentState(Map<String, String> values) {

        //TODO Go through generically and set the types
        String typeString = values.get("type");
        if (typeString != null) {
            type = Type.valueOf(typeString);
        }
        String minimumString = values.get("minimum");
        if (minimumString != null) {
            minimum = Float.parseFloat(minimumString);
        }
        String baselineString = values.get("baseline");
        if (typeString != null) {
            baseline = Float.parseFloat(baselineString);
        }
        String maximumString = values.get("maximum");
        if (typeString != null) {
            maximum = Float.parseFloat(maximumString);
        }
        String flipAxisString = values.get("flipAxis");
        if (typeString != null) {
            flipAxis = Boolean.parseBoolean(flipAxisString);
        }
        String drawBaselineString = values.get("drawBaseline");
        if (drawBaselineString != null) {
            drawBaseline = Boolean.parseBoolean(drawBaselineString);
        }
    }
}
