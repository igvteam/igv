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
package edu.mit.broad.prodinfo.datastrutures;

/**
 * Semi-open interval on the integer number line.
 * Turf covered runs from the start value inclusive, up to, but not including, the end value.
 *
 * @author tsharpe
 * @version $Revision: 33567 $
 */
public interface Interval
{
    // bit-wise definitions from which the other constants are composed
    static final int HAS_LESSER_PART = 1;
    static final int HAS_OVERLAPPING_PART = 2;
    static final int HAS_GREATER_PART = 4;

    static final int IS_ADJACENT_AND_EMPTY = 0;
    static final int IS_STRICTLY_LESS = HAS_LESSER_PART; // 1
    static final int IS_SUBSET = HAS_OVERLAPPING_PART; // 2
    static final int IS_LEFT_OVERHANGING_OVERLAPPER = HAS_LESSER_PART | HAS_OVERLAPPING_PART; // 3
    static final int IS_STRICTLY_GREATER = HAS_GREATER_PART; // 4
    // there is no value that equals 5, since that would imply overhanging on left and right without overlapping
    static final int IS_RIGHT_OVERHANGING_OVERLAPPER = HAS_GREATER_PART | HAS_OVERLAPPING_PART; // 6
    static final int IS_SUPERSET = HAS_LESSER_PART | HAS_OVERLAPPING_PART | HAS_GREATER_PART; // 7

    /**
     * Returns the starting point of the interval.
     * @return The start.
     */
    int getStart();

    /**
     * Returns the ending point of the interval.
     * The interval is not regarded as including this point.
     * @return The end.
     */
    int getEnd();

    /**
     * End - start.
     */
    int getLength();

    /**
     * Returns a constant that describes the relationship of this interval
     * to a specified interval with regard to position on the number line.
     * @param interval The interval to compare this one to.
     * @return One of the IS_* constants defined above.
     */
    int getRelationship( Interval interval );

    /**
     * Returns true if this interval ends where the specified interval starts,
     * or vice versa.
     * @param interval The interval to compare this one to.
     * @return True, if adjacent.
     */
    boolean isAdjacent( Interval interval );
}

