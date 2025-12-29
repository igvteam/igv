package org.broad.igv.feature;

import org.broad.igv.track.WindowFunction;

import java.awt.*;

/**
 * Representation of a feature from a UCSC "snp" file
 *
 * @author jrobinso
 *         Date: 11/5/13
 *         Time: 1:11 PM
 */
public class UCSCSnpFeature implements IGVFeature, htsjdk.tribble.Feature {

    String chr;
    int start;
    int end;
    float score;
    Strand strand;
    String name;
    String observed;
    String molType;
    String snpClass;
    String function;
    String submitters;
    String alleles;
    String alleleFreqs;

    public UCSCSnpFeature(String chr, int start, int end, String[] tokens) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.name = tokens[4];
        this.score = tokens[5].equals(".") ? 1000 : Float.parseFloat(tokens[5]);
        this.strand = tokens[6].equals("+") ? Strand.POSITIVE : (tokens[6].equals("-") ? Strand.NEGATIVE : Strand.NONE);
        this.observed = tokens[9];
        this.molType = tokens[10];
        this.snpClass = tokens[11];
        this.function = tokens[15];
        this.submitters = tokens[20];
        this.alleles = tokens[22];
        this.alleleFreqs = tokens[24];
    }

    @Override
    public String getChr() {
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

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public float getScore() {
        return this.score;
    }


    @Override
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        return getDescription();
    }


    public String getDescription() {

        StringBuffer desc = new StringBuffer();

        if (end - start > 1) {
            desc.append(chr + ":" + (start + 1) + "-" + end);
        } else {
            desc.append(chr + ":" + (start + 1));
        }

        desc.append("<br><b>Name:</b> " + name);
        desc.append("<br><b>Observed:</b> " + observed);
        desc.append("<br><b>Mol type:</b> " + molType);
        desc.append("<br><b>Class:</b> " + snpClass);
        desc.append("<br><b>Function:</b> " + function);
        desc.append("<br><b>Alleles:</b> " + alleles.replace(",", ", "));
        desc.append("<br><b>Allele freqs:</b> " + alleleFreqs.replace(",", ", "));
        desc.append("<br><b>Submitters:</b> " + submitters.replace(",", ", "));
        return desc.toString();
    }


    // Everything below has to be implemented for IGVFeature.   Sigh....

    @Override
    public Strand getStrand() {
        return strand;
    }

    @Override
    public Color getColor() {
        return Color.black;
    }

}
