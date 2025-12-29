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
public class Locus extends Range implements IGVNamedFeature {

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

    public void merge(Locus other) {
        if (this.getChr().equals(other.getChr())) {
            this.start = Math.min(this.start, other.start);
            this.end = Math.max(this.end, other.end);
        } else {
            throw new IllegalArgumentException("Cannot merge loci on different chromosomes");
        }
    }


    // Only accept full locus strings,  i.e. must contain : and -
    private void parseLocusString(String locusString) {
        int colonIndex = locusString.indexOf(":");
        int dashIndex = locusString.indexOf("-");
        if (colonIndex > 0) {
            chr = locusString.substring(0, colonIndex);
            String posString = locusString.substring(colonIndex).replace(":", "");
            setStartEnd(posString);
        }
    }

    private void setStartEnd(String posString) {
        String[] posTokens = posString.split("-");
        String startString = posTokens[0].replaceAll(",", "");
        start = Math.max(0, Integer.parseInt(startString) - 1);
        end = start + 1;
        if (posTokens.length >= 2) {
            try {
                String endString = posTokens[1].replaceAll(",", "");
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
        String startStr = NUMBER_FORMAT.format(start + 1);
        String endStr = NUMBER_FORMAT.format(end);
        String position = chr + ":" + startStr + "-" + endStr;
        return position;
    }

}
