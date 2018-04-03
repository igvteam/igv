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
