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

package org.broad.igv.util.index;


// Quick and dirty interval class
public class Interval implements Comparable {
    private final int low;
    private final int high;

    public Interval(int low, int high) {
        assert low <= high;
        this.low = low;
        this.high = high;
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


    public int getLow() {
        return this.low;
    }

    public int getHigh() {
        return this.high;
    }
}

