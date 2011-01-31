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
import org.broad.igv.sam.Alignment;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * A convenience class providing default implementation for many IGVFeature
 * methods.
 *
 * @author jrobinso
 */
public class BasicFeature extends AbstractFeature {

    private static Logger log = Logger.getLogger(BasicFeature.class);

    //A special, known feature type that represents a splice junction. A feature of this type has
    //slightly different popup text.
    //todo: refactor features and parsers to do this check in a cleaner way
    public static final String FEATURE_TYPE_SPLICEJUNCTION = "FEATURE_TYPE_SPLICEJUNCTION";

    protected List<Exon> exons;
    protected int level = 1;
    private String type;
    protected float score = Float.NaN;
    protected float confidence;
    String identifier;
    private int thickEnd;
    private int thickStart;

    public String[] getParentIds() {
        return parentIds;
    }

    public void setParentIds(String[] parentIds) {
        this.parentIds = parentIds;
    }

    String[] parentIds;
    String link;

    /**
     * Constructs ...
     */
    public BasicFeature() {
    }

    /**
     * Constructs ...
     *
     * @param chr
     * @param start
     * @param end
     */
    public BasicFeature(String chr, int start, int end) {

        this(chr, start, end, Strand.NONE);
        this.thickStart = start;
        this.thickEnd = end;
    }

    /**
     * Constructs ...
     *
     * @param chr
     * @param start
     * @param end
     * @param strand
     */
    public BasicFeature(String chr, int start, int end, Strand strand) {
        super(chr, start, end, strand);
        this.thickStart = start;
        this.thickEnd = end;

    }

    /**
     * Constructs ...
     *
     * @param feature
     */
    public BasicFeature(BasicFeature feature) {
        super(feature.getChr(), feature.getStart(), feature.getEnd(), feature.getStrand());
        super.setName(feature.getName());
        this.confidence = feature.confidence;
        this.color = feature.color;
        this.description = feature.description;
        this.exons = feature.exons;
        this.level = feature.level;
        this.score = feature.score;
        this.identifier = feature.identifier;
        this.type = feature.type;
        this.link = feature.link;
        this.thickStart = feature.thickStart;
        this.thickEnd = feature.thickEnd;
    }

    /**
     * Method description
     *
     * @return
     */
    public BasicFeature copy() {
        return new BasicFeature(this);
    }

    /**
     * @param identifier
     */
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }


    /**
     * Return a string for popup text, and related uses.  The default just
     * returns the feature name.  Its expected that this method will be
     * overriden in subclasses.
     *
     * @position -- 1 based coordinates
     */
    public String getValueString(double position, WindowFunction ignored) {
        StringBuffer valueString = new StringBuffer();

        // Get exon number, if over an exon
        int posZero = (int) position - 1;
        if (this.exons != null) {
            for (Exon exon : exons) {
                if (posZero >= exon.getStart() && posZero < exon.getEnd()) {
                    valueString.append(exon.getValueString(position, ignored) + "<br>");
                    return valueString.toString();
                }
            }
        }

        String name = getName();
        if (name != null) {
            valueString.append(name);
        }
        if ((identifier != null) && ((name == null) || !name.equals(identifier))) {
            valueString.append("<br>" + identifier);
        }
        valueString.append("<br>" + getLocusString());

        if (FEATURE_TYPE_SPLICEJUNCTION.equals(getType()) && getStrand() != null)
            valueString.append("<br>Strand: " + (getStrand().equals(Strand.POSITIVE) ? "+" : "-"));
        if (hasScore()) {
            //if this feature is a splice junction, then we're treating "Score" as "Depth".
            //todo: this is a hack.  Remove it after a refactor that allows it to be removed.
            String scoreDisplayName =
                    FEATURE_TYPE_SPLICEJUNCTION.equals(getType()) ? "Depth" : "Score";
            valueString.append("<br>" + scoreDisplayName + " = " + score);
        }
        if (description != null) {
            valueString.append("<br>" + description);
        }
        return valueString.toString();
    }


    /**
     * Method description
     *
     * @param score
     */
    public void setScore(float score) {
        this.score = score;
    }

    /**
     * Method description
     *
     * @return
     */
    @Override
    public float getScore() {
        return score;
    }


    @Override
    public boolean hasScore() {
        return Float.isNaN(score) ? false : true;
    }


    @Override
    public List<Exon> getExons() {
        return exons;
    }


    // TODO -- reimplement for efficienty

    @Override
    public Exon getExonAt(double position) {
        if (exons == null) return null;
        for (Exon exon : exons) {
            if (position >= exon.getStart() && position <= exon.getEnd()) {
                return exon;
            }
        }
        return null;
    }


    public void sortExons() {
        if (exons != null) {
            Collections.sort(exons, new Comparator<IGVFeature>() {

                public int compare(IGVFeature arg0, IGVFeature arg1) {
                    return arg0.getStart() - arg1.getStart();
                }
            });
        }
    }

    public void addExon(Exon region) {
        if (exons == null) {
            exons = new ArrayList();
        }
        setStart(Math.min(getStart(), region.getStart()));
        setEnd(Math.max(getEnd(), region.getEnd()));
        exons.add(region);
    }


    public int compareTo(Object o) {
        IGVFeature otherObj = (IGVFeature) o;
        return (int) (getStart() - otherObj.getStart());
    }

    /**
     * Method description
     *
     * @return
     */
    @Override
    public String getIdentifier() {
        return identifier;
    }

    /**
     * Method description
     *
     * @return
     */
    public int getExonCount() {
        return (exons == null) ? 0 : exons.size();
    }

    /**
     * Method description
     *
     * @return
     */
    public String getType() {
        return type;
    }

    /**
     * Method description
     *
     * @param type
     */
    public void setType(String type) {
        this.type = type;
    }

    public void setURL(String link) {
        this.link = link;
    }

    public String getURL() {
        return link;
    }


    public int getThickEnd() {
        return thickEnd;
    }

    public void setThickEnd(int thickEnd) {
        this.thickEnd = thickEnd;
    }

    public int getThickStart() {
        return thickStart;
    }

    public void setThickStart(int thickStart) {
        this.thickStart = thickStart;
    }
}
