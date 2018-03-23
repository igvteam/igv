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
package org.broad.igv.feature;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * @author jrobinso
 */
public class Locus extends Range implements NamedFeature {

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

    public Locus(String locusString) {
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
