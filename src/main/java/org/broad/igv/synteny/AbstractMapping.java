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

public abstract class AbstractMapping implements Mapping {
    private String name;
    protected String fromChr;
    protected int fromStart;
    protected int fromEnd;
    protected String toChr;
    protected int toStart;
    protected int toEnd;
    protected boolean direction;
    double scaleFactor;

    public void setParameters(String name, String fromChr, int fromStart, int fromEnd, String fromDir,
                              String toChr, int toStart, int toEnd, String toDir) {
        this.name = name;
        this.fromChr = fromChr;
        this.fromStart = fromStart;
        this.fromEnd = fromEnd;
        this.toChr = toChr;
        this.toStart = toStart;
        this.toEnd = toEnd;
        this.direction = fromDir.equals(toDir);
        this.scaleFactor = ((double) (toEnd - toStart) / (fromEnd - fromStart));
    }

    public String toString() {
        return this.name + " " + this.fromChr + ":" + this.fromStart + "-" + this.fromEnd + " -> " + this.toChr + ":" + this.toStart + "-" + this.toEnd;
    }

    public boolean containsFromPosition(int fromPosition) {
        return (fromPosition >= this.fromStart) && (fromPosition <= this.fromEnd);
    }

    public String getName() {
        return this.name;
    }

    public String getFromChr() {
        return this.fromChr;
    }

    public int getFromStart() {
        return this.fromStart;
    }

    public int getFromEnd() {
        return this.fromEnd;
    }

    public String getToChr() {
        return this.toChr;
    }

    public int getToStart() {
        return this.toStart;
    }

    public int getToEnd() {
        return this.toEnd;
    }

    public boolean getDirection() {
        return this.direction;
    }
}