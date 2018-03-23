package org.broad.igv.variant.New;

/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Broad Institute
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


import htsjdk.tribble.Feature;
import org.broad.igv.Globals;

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
