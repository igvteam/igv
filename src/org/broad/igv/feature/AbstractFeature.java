/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
import org.broad.igv.util.collections.MultiMap;
import org.broad.tribble.Feature;

import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
abstract public class AbstractFeature implements IGVFeature, org.broad.tribble.Feature {

    private static Logger log = Logger.getLogger(AbstractFeature.class);
    protected Strand strand = Strand.NONE;
    protected String chromosome;
    protected int start = -1;
    protected int end = -1;
    protected String type = "";
    protected Color color;
    protected String description;
    protected MultiMap<String, String> attributes;
    protected String name = "";

    /**
     * Constructs ...
     */
    public AbstractFeature() {
    }

    /**
     * Constructs ...
     *
     * @param chr
     * @param start
     * @param end
     * @param strand
     */
    public AbstractFeature(String chr, int start, int end, Strand strand) {
        this.chromosome = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }

    public String getIdentifier() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void setType(String type) {
        this.type = type;
    }

    public String getType() {
        return type;
    }

    public String getName() {
        return name;
    }

    public List<Exon> getExons() {
        return null;
    }

    public boolean hasExons() {
        return false;
    }

    public float getScore() {
        return Float.NaN;
    }

    public String getChr() {
        return chromosome;
    }

    /**
     * By default features are 1 bp wide
     *
     * @return
     */
    public int getEnd() {
        return end;
    }

    public int getLength() {
        return end - start;
    }

    public int getStart() {
        return start;
    }

    /**
     * Return true if the feature is completely contained within the bounds of this
     * feature. amd is on the same strand..
     * <p/>
     *
     * @param feature
     * @return
     */
    public boolean contains(IGVFeature feature) {
        if (feature == null) {
            return false;
        }
        if (!this.getChr().equals(feature.getChr()) ||
                this.getStrand() != feature.getStrand()) {
            return false;
        }
        if ((feature.getStart() >= this.getStart()) && (feature.getEnd() <= this.getEnd())) {
            return true;
        } else {
            return false;
        }
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public Strand getStrand() {
        return strand;
    }

    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    public boolean hasStrand() {
        return ((strand != null) && (strand != Strand.NONE));
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
    }

    public void setColor(String[] rgb, int nTokens) {
        try {
            if (nTokens < 3) {
                if (rgb[0].equals(".")) {
                    return;
                }
                color = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[0]),
                        Integer.parseInt(rgb[0]));
            } else {
                color = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]),
                        Integer.parseInt(rgb[2]));
            }
        } catch (NumberFormatException numberFormatException) {

        }
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getDescription() {
        return (description == null) ? getName() : description;
    }

    public MultiMap<String, String> getAttributes() {
        return attributes;
    }

    public void setAttributes(MultiMap<String, String> attributes) {
        this.attributes = attributes;
    }


    public void setAttribute(String key, String value) {
        if(attributes == null) {
            attributes = new MultiMap<String, String>();
        }
        attributes.put(key, value);
    }

    public boolean contains(double location) {
        return location >= getStart() && location < getEnd();
    }

    public boolean overlaps(Feature anotherFeature) {

        return end >= anotherFeature.getStart() && start <= anotherFeature.getEnd() &&
                chromosome.equals(anotherFeature.getChr());

    }


    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @param chromosome the chromosome to set
     */
    public void setChr(String chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * Get human readable locus string.
     * It is assumed that this.start and this.end are 0-based
     * and end-exclusive, and the returned string is 1-based
     * and end-inclusive. So basically we just add 1 to the start.
     *
     * @return
     */
    public String getLocusString() {
        return getChr() + ":" + (getStart() + 1) + "-" + getEnd();
    }


    protected String getAttributeString() {

        StringBuffer buf = new StringBuffer();
        buf.append("<br>");
        // 30 attributes is the maximum visible on a typical screen
        int max = 30;
        attributes.printHtml(buf, max);
        return buf.toString();

    }
}
