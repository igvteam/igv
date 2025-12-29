package org.igv.variant.New;

import htsjdk.tribble.Feature;
import org.igv.Globals;

/**
 * Parser for VCF files.
 */

public class Variant implements Feature {

    String chr;
    int pos;
    String names;
    String referenceBases;
    String alternateBases;
    int quality;
    String filter;
    String info;

    int start;
    int end;
    String[] alleles;


    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }



    public String getValueString(int genomicPosition) {

        StringBuffer buffer = new StringBuffer();

        buffer.append("Chr: " + this.chr + "<br>");
        buffer.append("Pos: " + (this.pos + 1));

        if (this.names != null) {
            buffer.append("Names: " + this.names + "<br>");
        }
        buffer.append("Ref: " + this.alleles[0] + "<br>");

        if (this.alleles.length > 1) {
            for (int i = 1; i < this.alleles.length; i++) {
                buffer.append(this.alleles[i]);
                if (i < this.alleles.length - 1) {
                    buffer.append(", ");
                }
                buffer.append("<br>");
            }
        }

        buffer.append("Quality: " + this.quality + "<br>");
        buffer.append("Filter: " + this.filter);

        if (this.info != null) {
            buffer.append("<hr>");
            String[] infoFields = Globals.semicolonPattern.split(this.info);

            for (String ifield : infoFields) {
                buffer.append(ifield + "<br>");

            }


        }

        return buffer.toString();
    }


}
