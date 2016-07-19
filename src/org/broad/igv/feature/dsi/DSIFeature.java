package org.broad.igv.feature.dsi;

import htsjdk.tribble.Feature;

/**
 * Created by jrobinson on 7/19/16.
 * <p>
 * Chr	 -	Number of chromosome
 * Position -	Start position
 * Base 	 -	Nucleobase cytosine or guanine
 * Total	 -	Total number of reads
 * Meth.	 -	Number of methylated reads
 * Unmeth.	 -	Number of unmethylation detection
 * Type	 -	CG/GC/other/unknown
 * F	 -	Fully methylated
 * P	 -	Hemimethylated on top (plus) strand
 * M	 -	Hemimethylated on bottom (minus) strand
 * U	 -	Unmethylated
 */
public class DSIFeature implements Feature {

    String chr;
    int position;
    char base;
    int total;
    int meth;
    int unmeth;
    String type;
    int f;
    int p;
    int m;
    int u;

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
        return position;
    }

    @Override
    public int getEnd() {
        return position + 1;
    }


    public String getValueString(double position, Object o) {

        StringBuffer buffer = new StringBuffer();
        buffer.append("Chromosome: " + chr +
                "<br>Position: " + (position + 1) +
                "<br>Type: " + type +
                "<br>Base: " + base +
                "<br>Total: " + total +
                "<br>Methylated: " + meth +
                "<br>Unmethylated: " + unmeth);


        if (f != Integer.MIN_VALUE) {
            buffer.append("<hr>");
            buffer.append("Fully methylated: " + f +
                    "<br>Hemimethylated on top (plus) strand: " + p +
                    "<br>Hemimethylated on bottom (minus) strand: " + m +
                    "<br>Unmethylated: " + u);
        }

        return buffer.toString();

    }
}
