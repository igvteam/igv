/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.feature;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.util.ColorTable;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.List;
import java.util.Map;

/**
 * Represents a mutation
 * // TODO -- refactor this to not implement "IGVFeature"
 *
 * @author jrobinso
 */
public class Mutation implements IGVFeature {

    private static Logger log = Logger.getLogger(Mutation.class);
    private static Map<String, Color> colors;

    private String sampleId;
    private String chr;
    private int start;
    private int end;
    private String type;
    private Color color;
    private Map<String, String> attributes;


    /**
     * Constructs ...
     *
     * @param runId
     * @param chromosome
     * @param start
     * @param end
     * @param type
     */
    public Mutation(String runId, String chromosome, int start, int end, String type) {
        this.sampleId = runId;
        this.chr = chromosome;
        this.start = start;
        this.end = end;
        this.type = type;
    }

    /**
     * Constructs ...
     *
     * @param mutation
     */
    public Mutation(Mutation mutation) {
        this.sampleId = mutation.sampleId;
        this.chr = mutation.chr;
        this.start = mutation.start;
        this.end = mutation.end;
        this.type = mutation.type;
        this.color = mutation.color;
    }


    public void setChr(String chr) {
        this.chr = chr;
    }

    public void setName(String name) {
        type = name;
    }


    /**
     * Method description
     *
     * @return
     */
    public Mutation copy() {
        return new Mutation(this);
    }

    /**
     * Method description
     *
     * @return
     */
    public String getSampleId() {
        return sampleId;
    }

    /**
     * Method description
     *
     * @return
     */
    public String getType() {
        return "mutation";
    }

    /**
     * Method description
     *
     * @return
     */
    public String getMutationType() {
        return type;
    }

    /**
     * Method description
     *
     * @return
     */
    public String getName() {
        return type.toString();
    }

    /**
     * Method description
     *
     * @return
     */
    public String getDescription() {
        return getName();
    }

    private String valueString;

    public String getValueString(double position, WindowFunction ignored) {

        if (valueString == null && attributes != null) {
            StringBuffer buf = new StringBuffer();
            buf.append("<html>Type: ");
            buf.append(type +"<br>");
            for (Map.Entry<String, String> entry : attributes.entrySet()) {
                buf.append(entry.getKey() + ": " + entry.getValue() + "<br>");
            }
            valueString = buf.toString();
        }
        return valueString;
    }

    public boolean hasScore() {
        return false;
    }

    public Strand getStrand() {
        return Strand.NONE;
    }

    public boolean overlaps(IGVFeature track) {
        return false;
    }

    public String getChr() {
        return chr;
    }

    public void setColor(Color color) {

        // Ignore
    }

    public Color getColor() {
        ColorTable colorTable = PreferenceManager.getInstance().getMutationColorScheme() ;
        Color c = colorTable.get(getMutationType());
        return c;
    }

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

    public float getScore() {
        return 0;
    }

    /**
     * Return true if the feature is completely contained within the bounds of this
     * featre.
     *
     * @param feature
     * @return
     */
    public boolean contains(IGVFeature feature) {

        if (feature == null || !this.getChr().equals(feature.getChr())) {
            return false;
        }
        if ((feature.getStart() >= this.getStart()) && (feature.getEnd() <= this.getEnd())) {
            return true;
        } else {
            return false;
        }
    }

    public boolean contains(double location) {
        return location >= start && location <= end;
    }

    public String getURL() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Exon getExonAt(double location) {
        return null;
    }

    public List<Exon> getExons() {
        return null;
    }

    public String getIdentifier() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public AminoAcidSequence getAminoAcidSequence(int exonIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getCdEnd() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getCdStart() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getLength() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public Map<String, String> getAttributes() {
        return null;
    }

    public void setAttributes(Map<String, String> attributes) {
        this.attributes = attributes;
    }
}
