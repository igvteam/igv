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
 * Date: Apr 5, 2010
 * Time: 4:00:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class WigItem {

    private static Logger log = Logger.getLogger(WigItem.class);

    private int itemIndex;         // wig section item index number
    private String chromosome;     // mChromosome name
    private int startBase;         // mStartBase base position for feature
    private int endBase;           // mEndBase base position for feature
    private float wigValue;        // wig value

    public WigItem(int itemIndex, String chromosome, int startBase, int endBase, float wigValue){

        this.itemIndex = itemIndex;
        this.chromosome = chromosome;
        this.startBase = startBase;
        this.endBase = endBase;
        this.wigValue = wigValue;
    }

    public int getItemNumber(){
        return itemIndex;
    }

    public String getChromosome() {
        return chromosome;
    }

    public int getStartBase() {
        return startBase;
    }

    public int getEndBase() {
        return endBase;
    }

    public float getWigValue() {
        return wigValue;
    }

     public void print(){
       log.debug("Wig item index " + itemIndex);
       log.debug("mChromosome name: " + chromosome);
       log.debug("mChromosome start base = " + startBase);
       log.debug("mChromosome end base = " + endBase);
       log.debug("Wig value: \n" + wigValue);
   }
}
