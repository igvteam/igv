package org.broad.igv.ucsc.bb;


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
public class BedData {

    // BBFile Table 12- BigBed data format
    public String chr;
    public int chromStart;         // starting base for item
    public int chromEnd;           // ending base for item
    public String restOfFields;    // string containing custom fields

    public BedData(String chr, int chromStart, int chromEnd, String restOfFields) {
        this.chr = chr;
        this.chromStart = chromStart;
        this.chromEnd = chromEnd;
        this.restOfFields = restOfFields;
    }

    public String toString() {
        return "" + chr + "\t" + chromStart + "\t" + chromEnd + "\t" + restOfFields;
    }


}
