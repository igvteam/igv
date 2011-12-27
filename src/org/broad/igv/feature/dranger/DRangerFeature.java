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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.dranger;

import org.broad.igv.feature.AbstractFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class DRangerFeature extends AbstractFeature {

    //   //num	chr1	str1	pos1	chr2	str2	pos2	tumreads	normreads	class	span	site1	site2	quality	score

    private int index;

    String chr2;
    int pos2;
    Strand str2;

    private int tumreads;
    private int normreads;
    private String featureClass;
    private int span;
    private String site1;
    private String site2;
    private int quality;
    private int score;


    public DRangerFeature() {

    }

    /**
     * Construct a new dRanger feature.  Arbitrarily sized at 10 bases so they are easier to see
     *
     * @param chr
     * @param pos1
     * @param strand1
     * @param chr2
     * @param pos2
     * @param str2
     */

    public DRangerFeature(String chr, int pos1, Strand strand1, String chr2, int pos2, Strand str2) {
        super(chr, pos1 - 5, pos1 + 5, strand1);
        this.chr2 = chr2;
        this.pos2 = pos2;
        this.str2 = str2;
    }

    public String getURL() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    // Score calibrated to = ~ 1,000 if 1/2 the reads support the rearrangment

    @Override
    public float getScore() {
        return score;

    }

    public LocusScore copy() {
        DRangerFeature newFeat = new DRangerFeature();
        newFeat.setChr(getChr());
        newFeat.setStart(getStart());
        newFeat.setEnd(getEnd());
        newFeat.chr2 = chr2;
        newFeat.pos2 = pos2;
        newFeat.str2 = str2;
        newFeat.tumreads = tumreads;
        newFeat.normreads = normreads;
        newFeat.featureClass = featureClass;
        newFeat.span = span;
        newFeat.site1 = site1;
        newFeat.site2 = site2;
        newFeat.quality = quality;
        newFeat.score = score;
        return newFeat;
    }

    public void setSite1(String site1) {
        this.site1 = site1;
    }


    public void setIndex(int index) {
        this.index = index;
    }

    public void setTumreads(int tumreads) {
        this.tumreads = tumreads;
    }

    public void setNormreads(int normreads) {
        this.normreads = normreads;
    }

    public void setFeatureClass(String featureClass) {
        this.featureClass = featureClass;
    }


    public void setSpan(int span) {
        this.span = span;
        //setEnd(getStart() + span);
    }

    public void setQuality(int quality) {
        this.quality = quality;
    }

    public void setScore(int score) {
        this.score = score;
    }


    Color defaultColor = Color.blue;
    static Map<String, Color> colorMap;
    static {
        colorMap = new HashMap();
        colorMap.put("deletion", Color.red);
        colorMap.put("inter_chr", Color.blue);
        colorMap.put("inversion", Color.green);
        colorMap.put("long_range", Color.MAGENTA);
        colorMap.put("tandem_dup", Color.cyan);
    }

    @Override
    public Color getColor() {
        Color baseColor = defaultColor;
        if(featureClass != null  && colorMap.containsKey(featureClass)) {
            baseColor = colorMap.get(featureClass);
        }
        float alpha = Math.min(1.0f, Math.max(0.15f, .15f + (.85f / 3) * score));
        return ColorUtilities.getCompositeColor(baseColor, alpha);
    }


    /*
       int index;
   int tumorReads;
   int normalReads;
   private String featureClass;
   private int span;
   private String site = "";
   int quality;
   int score;
    */
    //num	chr1	str1	pos1	chr2	str2	pos2	tumreads	normreads	class	span

    // site1	site2	quality	score

    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer buf = new StringBuffer();
        if (index > 0) {
            buf.append("Index:   " + index + "<br>");
        }
        if (site1 != null) {
            buf.append("Site1:    " + site1 + "<br>");
        }
        if (site2 != null) {
            buf.append("Site2:    " + site2 + "<br>");
        }
        if (tumreads > 0) {
            buf.append("T reads: " + tumreads + "<br>");
        }
        if (normreads > 0) {
            buf.append("N reads: " + normreads + "<br>");
        }
        if (featureClass != null) {
            buf.append("Class: " + featureClass + "<br>");
        }
        if (span > 0) {
            buf.append("Span:    " + span + "<br>");
        }
        if (quality > 0) {
            buf.append("Quality:   " + quality + "<br>");
        }
        if (score > 0) {
            buf.append("Score:   " + score);
        }

        return buf.toString();

    }

    public void setSite2(String site2) {
        this.site2 = site2;
    }
}
