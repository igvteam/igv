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

package org.broad.igv.util.index;


// Quick and dirty interval class
public class Interval implements Comparable {

     final int low;
     final int high;
     final long value;

    public Interval(int low, int high, long value) {
        assert low <= high;
        this.low = low;
        this.high = high;
        this.value = value;
    }


    public boolean equals(Object other) {
        if (other == null) {
            return false;
        }
        if (this == other)
            return true;
        if (this.getClass().equals(other.getClass())) {
            Interval otherInterval = (Interval) other;
            return (this.low == otherInterval.low &&
                    this.high == otherInterval.high);
        }
        return false;
    }


    public int hashCode() {
        return low;
    }


    public int compareTo(Object o) {
        Interval other = (Interval) o;
        if (this.low < other.low)
            return -1;
        if (this.low > other.low)
            return 1;

        if (this.high < other.high)
            return -1;
        if (this.high > other.high)
            return 1;

        return 0;
    }

    public String toString() {
        return "Interval[" + this.low + ", " + this.high + "]";
    }


    /**
     * Returns true if this interval overlaps the other.
     */
    public boolean overlaps(Interval other) {
        return (this.low <= other.high &&
                other.low <= this.high);
    }

    public boolean overlaps(int low, int high) {
        return (this.low <= high && low <= this.high);

    }


    public int getLow() {
        return this.low;
    }

    public int getHigh() {
        return this.high;
    }

    public long getValue() {
        return value;
    }
}

