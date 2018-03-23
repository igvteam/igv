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
