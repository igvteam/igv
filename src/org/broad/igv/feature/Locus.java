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

import org.broad.igv.data.Interval;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * @author jrobinso
 */
public class Locus extends Interval implements NamedFeature {

    private static NumberFormat NUMBER_FORMAT = NumberFormat.getInstance(Locale.US);

    public Locus(String chr, int start, int end) {
        super(chr, start, end);
    }

    /**
     * Generate Locus from locusString of form [chr]:[start]-[end]
     * Returns null if invalid string
     * @param locusString
     * @return
     */
    public static Locus fromString(String locusString) {
        Locus l = new Locus(locusString);
        return l.isValid() ? l : null;
    }

    private Locus(String locusString) {
        this(null, -1, -1);
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

}
