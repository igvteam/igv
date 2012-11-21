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
package org.broad.igv.feature;

import org.broad.tribble.Feature;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * @author jrobinso
 */
public class Locus implements NamedFeature {

    protected String chr = null;
    protected int start = -1;
    protected int end = -1;

    private static NumberFormat NUMBER_FORMAT = NumberFormat.getInstance(Locale.US);

    public Locus(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

    public static Locus fromString(String locusString) {
        Locus l = new Locus(locusString);
        return l.isValid() ? l : null;
    }

    // TODO -- really the factory method above should be used and this constructor made private or removed
    public Locus(String locusString) {
        parseLocusString(locusString);

    }

    public boolean isValid() {
        return getChr() != null && getStart() >= 0 && getEnd() >= getStart();
    }

    // Only accept full locus strings,  i.e. must contain : and -
    private void parseLocusString(String locusString) {
        int colonIndex = locusString.indexOf(":");
        int dashIndex = locusString.indexOf("-");
        if (colonIndex > 0 && dashIndex > 0) {
            chr = locusString.substring(0, colonIndex);
            String posString = locusString.substring(colonIndex).replace(":", "");
            setStartEnd(posString);
        }
    }

    private void setStartEnd(String posString) {
        String[] posTokens = posString.split("-");
        if (posTokens.length >= 2) {
            try {
                String startString = posTokens[0].replaceAll(",", "");
                String endString = posTokens[1].replaceAll(",", "");
                start = Math.max(0, Integer.parseInt(startString));
                end = Integer.parseInt(endString);
            } catch (NumberFormatException numberFormatException) {
                // This can happen, indicates an invalid or incomplete locus string
                // For now nothing to do here, 

            }
        }
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public String toString() {
        return chr + ":" + start + "-" + end;
    }

    @Override
    public String getName() {
        return toString();
    }

    public static String getFormattedLocusString(String chr, int start, int end) {
        String startStr = NUMBER_FORMAT.format(start);
        String endStr = NUMBER_FORMAT.format(end);
        String position = chr + ":" + startStr + "-" + endStr;

        return position;
    }

    /**
     * Weak version of contains, does not check genome
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public boolean contains(String chr, int start, int end) {
        return this.chr.equals(chr) && this.start <= start && this.end >= end;
    }

    public boolean overlaps(String chr, int start, int end) {
        return this.chr.equals(chr) && this.start <= end && this.end >= start;
    }


}
