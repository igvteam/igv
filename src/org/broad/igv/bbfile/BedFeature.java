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
 * Date: Mar 18, 2010
 * Time: 3:36:10 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for BigBed features.
*
*   Note: required BigBed data items are:
*       mChromosome (name)
*       mChromosome mStartBase (starting base)
*       mChromosome mEndBase (ending base)
*       plus String "rest of fields" for custom fileds
*
*   Custom fields can follow any of the predefined fields which normally
*   follow the three required fields. Predefined fileds must be maintained
*   up to the point of customization.  (See BBFile Table A)
*
*   The predefined fields are:
*       name - name of feature
*       score - value betwenn 0 and 1000 defining viewing darkness
*       strand - "+" or "-", or "." for unknown
*       thickstart - base where item thickens, used for CDS mStartBase of genes
*       thickEnd - base where thick item ends
*       itemRGB - comma seperated R,G,B valuse from 0 to 255
*       blockCount - number of multi-part blocks; number of exons for genes
*       blockSizes - blockCount comma seperated list of blocks
*       blockStarts - blockCount comma seperated mStartBase locations (relative to mChromosome mStartBase)
*
*       Custom field dimensions are defined by the following fileds in BBFile Tab;e C:
 *          field count - number of fields in Bed format
 *          defined field count - number of fields that are of predefied type as shown above
*
*   Custom fields:
*       restOfFields (String contains the predefined and custom fields)
*
*   The custom fields are described by  .as dictionary terms which are
*   provided by the autoSQL section of the BigBed file. (See BBFile Table B example)
*
* */
public class BedFeature {

    private static Logger log = Logger.getLogger(BedFeature.class);

    private int itemIndex;     // data record index

    // BBFile Table I - BigBed data format
    private String chromosome;      // mChromosome/contig name
    private int startBase;         // starting base for item
    private int endBase;           // ending base for item
    private String[] restOfFields;    // string containing custom fields

    public BedFeature(int itemIndex, String chromosome, int startBase, int endBase, String restOfFieldsString){

       this.itemIndex = itemIndex;
       this.chromosome =  chromosome;
       this.startBase =  startBase;
       this.endBase = endBase;
       restOfFields = restOfFieldsString == null ? null : restOfFieldsString.split("\t");
   }

   // returns the data record index
   public int getItemIndex() {
       return itemIndex;
   }

   // returns the mChromosome ID (0, 1, etc.)
   public String getChromosome() {
       return chromosome;
   }

   // returns the mChromosome mStartBase base position
   public int getStartBase(){
       return startBase;
   }

   // returns the mChromosome mEndBase base position
   public int getEndBase() {
       return endBase;
   }

    public String[] getRestOfFields(){
        return restOfFields;
    }

   public void print(){

       System.out.println("BigBed feature item " + itemIndex);
       System.out.println("mChromosome name: " + chromosome);
       System.out.println("mChromosome start base= " + startBase);
       log.debug("mChromosome end base = " + endBase);
       System.out.println("Rest of fields: \n" + restOfFields);
   }
}
