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

package org.broad.igv.synteny;

/**
 * User: jrobinso
 * Date: Feb 19, 2010
 */ // region R:chr1:chr2:D1         chr1 47243     59894 +   chr2 111274749 111285190 +
// anchor A:chr1:3928:chr4:34422 chr1 17429161 17429302 + chr4 140098907 140099048 - 96.7
public class SyntenyMapping {

    private String name;
    private String fromChr;
    private int fromStart;
    private int fromEnd;
    private String toChr;
    private int toStart;
    private int toEnd;
    private boolean direction;
    double scaleFactor;


    public SyntenyMapping(String name, String fromChr, int fromStart, int fromEnd, String fromDir,
                          String toChr, int toStart, int toEnd, String toDir) {
        this.name = name;
        this.fromChr = fromChr;
        this.fromStart = fromStart;
        this.fromEnd = fromEnd;
        this.toChr = toChr;
        this.toStart = toStart;
        this.toEnd = toEnd;
        this.direction = fromDir.equals(toDir);
        this.scaleFactor = ((double) (toEnd - toStart)) / (fromEnd - fromStart);

    }


    public double mapPosition(int fromPosition) {
        return (direction == true)
                ? toStart + scaleFactor * (fromPosition - fromStart) :
                toEnd - scaleFactor * (fromPosition - fromStart);

    }

    public boolean containsFromPosition(int fromPosition) {
        return fromPosition >= fromStart && fromPosition < fromEnd;
    }

    public boolean containsFromInterval(int start, int end) {
        return end >= fromStart && start < fromEnd;
    }


    public String getName() {
        return name;
    }


    public String getFromChr() {
        return fromChr;
    }


    public int getFromStart() {
        return fromStart;
    }


    public int getFromEnd() {
        return fromEnd;
    }


    public String getToChr() {
        return toChr;
    }


    public int getToStart() {
        return toStart;
    }


    public int getToEnd() {
        return toEnd;
    }


    public boolean getDirection() {
        return direction;
    }

}
