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

package org.broad.igv.feature.tribble;

import org.broad.igv.feature.Exon;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.MultiMap;

import java.awt.*;
import java.util.*;

/**
 * @author jrobinso
 */
public class VCFFeature implements IGVFeature, org.broad.tribble.Feature {

    private String chr;
    private int start;
    private String id = "";
    private String name;
    private String ref;
    private String alt;
    private float quality;
    private String filter;
    private String info;

    public VCFFeature(String chr, int start, String id, String ref, String alt, float quality, String filter, String info) {
        this.chr =  chr;
        this.start = start;
        this.id = id;
        this.ref = ref;
        this.alt = alt;
        this.name = "";
        this.quality = quality;
        this.filter = filter;
        this.info = info.replace(";", "<br>");
    }

    public float getScore() {
        return getQuality();
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer buf = new StringBuffer();
        if (id != null && id.length() > 0) {
            buf.append(id + "<br>");
        }
        buf.append(ref + "->" + alt +  "<br>");
        buf.append("QUAL = " + quality + "<br>");
        buf.append("FILTER = " + filter + "<br>");
        buf.append("--------------<br>");
        buf.append(info);
        return buf.toString();

    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return start + 1;
    }

    public String getIdentifier() {
        return id;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getName() {
        return name;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getDescription() {
        return name;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean hasScore() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public LocusScore copy() {
        throw new java.lang.UnsupportedOperationException("Copy not supported");
    }

    public Strand getStrand() {
        return Strand.NONE;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean hasStrand() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setName(String name) {
        throw new java.lang.UnsupportedOperationException("Copy not supported");
    }

    public Color getColor() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean contains(IGVFeature feature) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public java.util.List<Exon> getExons() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getLength() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public MultiMap<String, String> getAttributes() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean contains(double location) {
        return location >= start && location < start + 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getURL() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Exon getExonAt(double location) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setEnd(int end) {
        throw new java.lang.UnsupportedOperationException("setEnd not supported");

    }

    public String getType() {
        return "VCF";  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        throw new java.lang.UnsupportedOperationException("setEnd not supported");
    }

    public float getQuality() {
        return quality;
    }


    public String getFilter() {
        return filter;
    }

    public String getInfo() {
        return info;
    }


    // TODO -- move this up
    public String getLocusString() {
        return getChr() + ":" + getStart() + ":" + getEnd();
    }

}
