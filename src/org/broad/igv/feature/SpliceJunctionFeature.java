/**
 * Copyright (c) 2010-2011 by Fred Hutchinson Cancer Research Center.  All Rights Reserved.

 * This software is licensed under the terms of the GNU Lesser General
 * Public License (LGPL), Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.

 * THE SOFTWARE IS PROVIDED "AS IS." FRED HUTCHINSON CANCER RESEARCH CENTER MAKES NO
 * REPRESENTATIONS OR WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED,
 * INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS,
 * WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL FRED HUTCHINSON CANCER RESEARCH
 * CENTER OR ITS TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR
 * ANY DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR
 * CONSEQUENTIAL DAMAGES, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS,
 * REGARDLESS OF  WHETHER FRED HUTCHINSON CANCER RESEARCH CENTER SHALL BE ADVISED,
 * SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

package org.broad.igv.feature;

import org.broad.igv.track.WindowFunction;


/**
 * A feature class for splice junctions, with depth information for flanking regions if available.
 * <p/>
 * The interpretation of 'start' and 'end' for these features is the start and end of the flanking regions.
 * This is a bit counterintuitive, but it's required for popup text to display correctly in the feature track
 *
 * @author dhmay
 */
public class SpliceJunctionFeature extends BasicFeature {
    //depth of coverage for the splice junction
    protected int junctionDepth = 0;

    //start and end locations of the junction
    protected int junctionStart = 0;
    protected int junctionEnd = 0;

    int[] startFlankingRegionDepthArray, endFlankingRegionDepthArray;

    public SpliceJunctionFeature(String chr, int start, int end) {
        super(chr, start, end);
        junctionStart = start;
        junctionEnd = end;
    }

    public SpliceJunctionFeature(String chr, int start, int end, Strand strand) {
        super(chr, start, end, strand);
        junctionStart = start;
        junctionEnd = end;
    }

    public SpliceJunctionFeature(String chr, int start, int end, Strand strand,
                                 int junctionStart, int junctionEnd) {
        this(chr, start, end, strand);
        this.junctionStart = junctionStart;
        this.junctionEnd = junctionEnd;
    }

    /**
     * Does this splice junction feature represent the same splice junction as another feature? This is used for
     * splice junction feature selection.
     *
     * @param otherFeature
     * @return
     */
    public boolean isSameJunction(SpliceJunctionFeature otherFeature) {
        if (otherFeature.getJunctionStart() == getJunctionStart() &&
                otherFeature.getJunctionEnd() == getJunctionEnd())
            return true;
        return false;
    }

    /**
     * "Adds a read" to the junction. We don't actually track all the reads associated with a junction
     * (though that might be interesting, in the future), but we track depth of coverage information. By
     * adding all the reads one by one, we can build up the depth of coverage map for the flanking regions
     * and the junction itself.
     *
     * @param readStart
     * @param readEnd
     */
    public void addRead(int readStart, int readEnd) {
        junctionDepth++;

        int newStartFlankingRegionSize = junctionStart - readStart;
        if (readStart < start) {
            int[] newStartFlankArray = new int[newStartFlankingRegionSize];
            if (startFlankingRegionDepthArray != null) {
                int offset = newStartFlankingRegionSize - getStartFlankingRegionLength();
                System.arraycopy(startFlankingRegionDepthArray, 0, newStartFlankArray,
                        offset, getStartFlankingRegionLength());
            }
            startFlankingRegionDepthArray = newStartFlankArray;
            start = readStart;
        }
        for (int i = getStartFlankingRegionLength() - newStartFlankingRegionSize;
             i < getStartFlankingRegionLength(); i++)
            startFlankingRegionDepthArray[i] = startFlankingRegionDepthArray[i] + 1;

        int newEndFlankingRegionSize = readEnd - junctionEnd;
        if (readEnd > end) {
            int[] newEndFlankArray = new int[newEndFlankingRegionSize];
            if (endFlankingRegionDepthArray != null) {
                System.arraycopy(endFlankingRegionDepthArray, 0, newEndFlankArray,
                        0, getEndFlankingRegionLength());
            }
            endFlankingRegionDepthArray = newEndFlankArray;
            end = readEnd;
        }
        for (int i = 0; i < newEndFlankingRegionSize; i++)
            endFlankingRegionDepthArray[i] = endFlankingRegionDepthArray[i] + 1;
    }

    public int getJunctionDepth() {
        return junctionDepth;
    }


    public void setJunctionDepth(int junctionDepth) {
        this.junctionDepth = junctionDepth;
    }

    public int getJunctionStart() {
        return junctionStart;
    }

    public void setJunctionStart(int junctionStart) {
        this.junctionStart = junctionStart;
    }

    public int getJunctionEnd() {
        return junctionEnd;
    }

    public void setJunctionEnd(int junctionEnd) {
        this.junctionEnd = junctionEnd;
    }

    public int getStartFlankingRegionLength() {
        return junctionStart - start;
    }

    public int getEndFlankingRegionLength() {
        return end - junctionEnd;
    }

    public int[] getStartFlankingRegionDepthArray() {
        return startFlankingRegionDepthArray;
    }

    public int[] getEndFlankingRegionDepthArray() {
        return endFlankingRegionDepthArray;
    }

    public boolean hasFlankingRegionDepthArrays() {
        return ((startFlankingRegionDepthArray != null) && (endFlankingRegionDepthArray != null));
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
        String name = getName();
        if (name != null) {
            valueString.append(name);
        }
        if ((identifier != null) && ((name == null) || !name.equals(identifier))) {
            valueString.append("<br>" + identifier);
        }

        valueString.append("<br>");
        valueString.append(chromosome + ":" + junctionStart + "-" + junctionEnd);

        valueString.append("<br>Strand: " + (getStrand().equals(Strand.POSITIVE) ? "+" : "-"));
        valueString.append("<br>Depth = " + junctionDepth + ", Flanking Widths: (" +
                this.getStartFlankingRegionLength() + "," +
                this.getEndFlankingRegionLength() + ")");
        if (hasFlankingRegionDepthArrays()) {
            if (position >= start && position < junctionStart) {
                int index = (int) position - start;
                if (index < startFlankingRegionDepthArray.length)
                    valueString.append("<br>Start Flanking, Depth = " +
                            startFlankingRegionDepthArray[index]);
            } else if (position > junctionEnd && position <= end) {
                int index = (int) position - junctionEnd;
                if (index < endFlankingRegionDepthArray.length)
                    valueString.append("<br>End Flanking, Depth = " +
                            endFlankingRegionDepthArray[index]);
            }
        }
        if (description != null) {
            valueString.append("<br>" + description);
        }
        return valueString.toString();
    }

}
