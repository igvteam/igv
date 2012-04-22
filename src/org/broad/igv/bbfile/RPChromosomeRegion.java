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

package org.broad.igv.bbfile;

import org.apache.log4j.Logger;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 22, 2010
 * Time: 3:37:38 PM
 * private int mStartChromID;  // starting chromosome in item
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for R+ Tree bounding rectangle regions
* */
public class RPChromosomeRegion {

    private static Logger log = Logger.getLogger(RPChromosomeRegion.class);

    private int startChromID;  // starting mChromosome in item
    private int startBase;     // starting base pair in item
    private int endChromID;    // ending mChromosome in item
    private int endBase;       // ending base pair in item

    /*
    *   Construct region from a specification.
    * */

    public RPChromosomeRegion(int startChromID, int startBase,
                              int endChromID, int endBase) {

        this.startChromID = startChromID;
        this.startBase = startBase;
        this.endChromID = endChromID;
        this.endBase = endBase;
    }

    /*
   *   Construct region from an existing region.
   * */

    public RPChromosomeRegion(RPChromosomeRegion region) {

        startChromID = region.startChromID;
        startBase = region.startBase;
        endChromID = region.endChromID;
        endBase = region.endBase;
    }

    /*  
    *   Null region constructor for setting region members.
    **/

    public RPChromosomeRegion() {
        // members auto-inited
    }

    public int getStartChromID() {
        return startChromID;
    }

    public int getStartBase() {
        return startBase;
    }


    public int getEndChromID() {
        return endChromID;
    }


    public int getEndBase() {
        return endBase;
    }


    public void print() {

        log.debug("Chromosome bounds:");
        log.debug("StartChromID = " + startChromID);
        log.debug("StartBase = " + startBase);
        log.debug("EndChromID = " + endChromID);
        log.debug("EndBase = " + endBase);
    }

    /**
     * Comparator for mChromosome bounds is used to find relevant intervals and
     * rank placement of node items. Returned value indicates relative
     * positioning to supplied chromosome test region , and expands on normal
     * comparator by indicating partial overlap in the extremes.
     * <p/>
     * Returns:
     * - 2 indicates that this region is completely disjoint below the test region
     * -1 indicates this region intersects the test region from below
     * 0 indicates that this region is inclusive to the test region
     * 1 indicates this region intersects the test region from above
     * 2 indicates that this region is completely disjoint above the test region
     * <p/>
     * Note: additional tests can be applied to determine intersection from above
     * or below the test region and disjoint above or below the test region cases.
     */

    public int compareRegions(RPChromosomeRegion testRegion) {

        return compareRegions(testRegion.startChromID, testRegion.startBase, testRegion.endChromID, testRegion.endBase);
    }

    public int compareRegions(int testRegionStartChromID, int testRegionStartBase, int testRegionEndChromID, int testRegionEndBase) {

        // test if this region is contained by (i.e. subset of) testRegion region
        if (containedIn(testRegionStartChromID, testRegionStartBase, testRegionEndChromID, testRegionEndBase))
            return 0;

            // test if  testRegion region is disjoint from above or below
        else if (disjointBelow(testRegionStartChromID, testRegionStartBase))
            return -2;
        else if (disjointAbove(testRegionEndChromID, testRegionEndBase))
            return 2;

            // Otherwise this region must intersect
        else if (this.intersectsBelow(testRegionStartChromID, testRegionStartBase))
            return -1;
        else if (this.intersectsAbove(testRegionEndChromID, testRegionEndBase))
            return 1;

        // unexpected condition is unknown
        return 3;
    }


    /**
     * Method checks if test region contains this region;
     * (i.e this region is subset oftest region).
     * <p/>
     * Parameters:
     * testRegion - chromosome selection region
     * <p/>
     * Returns:
     * This region is contained in the test region: true or false
     */


    private boolean containedIn(int testRegionStartChromID, int testRegionStartBase, int testRegionEndChromID, int testRegionEndBase) {


        if (startChromID > testRegionStartChromID ||
                (startChromID == testRegionStartChromID && startBase >= testRegionStartBase)) {
            if (endChromID < testRegionEndChromID ||
                    (endChromID == testRegionEndChromID && endBase <= testRegionEndBase))
                return true;
            else
                return false;
        } else
            return false;
    }


    /**
     * Method checks if this region intersects test region from below
     * <p/>
     * Note: To be true, this region must have some part outside the test region
     * <p/>
     * Parameters:
     * testRegion - chromosome selection region
     * <p/>
     * Returns:
     * This region intersects the test region from below: true or false
     */
    private boolean intersectsBelow(int testRegionStartChromID, int testRegionStartBase) {

        // Only need to test if some part of this region is below and some within test region.
        if (startChromID < testRegionStartChromID ||
                (startChromID == testRegionStartChromID && startBase < testRegionStartBase)) {
            if (endChromID > testRegionStartChromID ||
                    (endChromID == testRegionStartChromID && endBase > testRegionStartBase))
                return true;
            else
                return false;
        } else
            return false;
    }

    /**
     * Method checks if this region intersects test region from above.
     * <p/>
     * Note: To be true, this region must have some part outside the test region
     * <p/>
     * Parameters:
     * testRegion - chromosome selection region
     * <p/>
     * Returns:
     * This region intersects the test region from above: true or false
     */
    private boolean intersectsAbove(int testRegionEndChromID, int testRegionEndBase) {

        // Only need to test if some part of this region is above and some within test region.
        if (endChromID > testRegionEndChromID ||
                (endChromID == testRegionEndChromID && endBase > testRegionEndBase)) {
            if (startChromID < testRegionEndChromID ||
                    (startChromID == testRegionEndChromID && startBase < testRegionEndBase))
                return true;
            else
                return false;
        } else
            return false;
    }

    /**
     * Method checks if this region is completely below test region.
     * <p/>
     * Parameters:
     * testRegion - chromosome selection region
     * <p/>
     * Returns:
     * This region is disjoint below the test region: true or false
     */
    private boolean disjointBelow(int testRegionStartChromID, int testRegionStartBase) {

        if (endChromID < testRegionStartChromID ||
                endChromID == testRegionStartChromID && endBase <= testRegionStartBase)
            return true;
        else
            return false;
    }

    /**
     * Method checks if this region region is completely above test region.
     * <p/>
     * Parameters:
     * testRegion - chromosome selection region
     * <p/>
     * Returns:
     * This region is disjoint above the test region: true or false
     */
    private boolean disjointAbove(int testRegionEndChromID, int testRegionEndBase) {

        if (startChromID > testRegionEndChromID ||
                startChromID == testRegionEndChromID && startBase >= testRegionEndBase)
            return true;
        else
            return false;
    }

    /**
     * Method computes the extremes between this region and the test region
     * <p/>
     * Parameters:
     * testRegion - chromosome region to compare against this region
     * <p/>
     * Returns:
     * new chromosome region of extremes
     */

    public RPChromosomeRegion getExtremes(RPChromosomeRegion testRegion) {
        RPChromosomeRegion newRegion = new RPChromosomeRegion(this);

        // update node bounds
        if (testRegion.startChromID < newRegion.startChromID ||
                (testRegion.startChromID == newRegion.startChromID &&
                        testRegion.startBase < newRegion.startBase)) {
            newRegion.startChromID = testRegion.startChromID;
            newRegion.startBase = testRegion.startBase;
        }

        if (testRegion.endChromID > newRegion.endChromID ||
                (testRegion.endChromID == newRegion.endChromID &&
                        testRegion.endBase > newRegion.endBase)) {
            newRegion.endChromID = testRegion.endChromID;
            newRegion.endBase = testRegion.endBase;
        }

        return newRegion;
    }

}
