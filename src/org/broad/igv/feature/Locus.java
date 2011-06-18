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
package org.broad.igv.feature;

/**
 * @author jrobinso
 */
public class Locus {

    private String chr = null;
    protected int start = -1;
    protected int end = -1;

    public Locus(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

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

    protected boolean overlaps(String chr, int start, int end) {
        return this.chr.equals(chr) && this.start <= end && this.end >= start;
    }
}
