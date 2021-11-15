package org.broad.igv.variant.vcf;

import org.broad.igv.variant.Allele;
import org.broad.igv.variant.Genotype;
import org.broad.igv.variant.Variant;

import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * Represents the mate of a structural variant, defined by CHR2 and END attributes.
 */
public class MateVariant implements Variant {

    private String chr;
    private int position;
    public Variant mate;

    public MateVariant(String chr, int position, Variant mate) {
        this.chr = chr;
        this.position = position;
        this.mate = mate;
    }

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public int getStart() {
        return position - 1;
    }

    @Override
    public int getEnd() {
        return position;
    }


    @Override
    public String getID() {
        return mate.getID();
    }

    @Override
    public boolean hasLog10PError() {
        return mate.hasLog10PError();
    }

    @Override
    public String getType() {
        return mate.getType();
    }

    @Override
    public boolean isFiltered() {
        return mate.isFiltered();
    }

    @Override
    public double getPhredScaledQual() {
        return mate.getPhredScaledQual();
    }

    @Override
    public Map<String, Object> getAttributes() {
        return mate.getAttributes();
    }

    @Override
    public String getAttributeAsString(String key) {
        return mate.getAttributeAsString(key);
    }

    @Override
    public String getReference() {
        return mate.getReference();
    }

    @Override
    public List<Allele> getAlternateAlleles() {
        return mate.getAlternateAlleles();
    }

    @Override
    public double[] getAlleleFreqs() {
        return mate.getAlleleFreqs();
    }

    @Override
    public Collection<String> getSampleNames() {
        return mate.getSampleNames();
    }

    @Override
    public Genotype getGenotype(String sample) {
        return mate.getGenotype(sample);
    }

    @Override
    public Collection<String> getFilters() {
        return mate.getFilters();
    }

    @Override
    public double getMethlationRate() {
        return mate.getMethlationRate();
    }

    @Override
    public double getCoveredSampleFraction() {
        return mate.getCoveredSampleFraction();
    }

    @Override
    public String getPositionString() {
        return mate.getPositionString();
    }

    @Override
    public int[] getAlleleCounts() {
        return mate.getAlleleCounts();
    }

    @Override
    public double getAlternateAlleleFrequency() {
        return mate.getAlternateAlleleFrequency();
    }

    @Override
    public int getTotalAlleleCount() {
        return mate.getTotalAlleleCount();
    }

    @Override
    public double getAlleleFraction() {
        return mate.getAlleleFraction();
    }

}
