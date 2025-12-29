package org.igv.variant.vcf;

import org.igv.logging.*;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.variant.Allele;
import org.igv.variant.Genotype;
import org.igv.variant.Variant;
import org.igv.variant.VariantTrack;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.text.DecimalFormat;
import java.util.*;

/**
 * @author Jim Robinson, jacob
 * @date Aug 1, 2011
 */
public class VCFVariant implements Variant {

    private static Logger log = LogManager.getLogger(Variant.class);

    static final DecimalFormat numFormat = new DecimalFormat("###,###,###");

    VariantContext variantContext;
    List<Allele> alternateAlleles;
    // private ZygosityCount zygosityCount;

    String chr;
    private double[] alleleFreqs;
    private int[] alleleCounts;
    private double methylationRate = Double.NaN;  // <= signals unknown / not applicable
    private double coveredSampleFraction = Double.NaN;

    Map<String, VCFGenotype> genotypeMap;

    private int start = -1;
    private int totalAlleleCount = 0;

    private boolean nonRef;

    public VCFVariant(VariantContext variantContext, String chr) {
        this.variantContext = variantContext;
        this.chr = chr;
        init();
    }

    private void init() {

        // Copy the genotype map.  Calls to variantContext.getGenotype() are expensive
        genotypeMap = new HashMap<String, VCFGenotype>();
        for (String sample : getSampleNames()) {
            htsjdk.variant.variantcontext.Genotype genotype = variantContext.getGenotype(sample);
            VCFGenotype vcfGenotype = genotype == null ? null : new VCFGenotype(genotype);
            genotypeMap.put(sample, vcfGenotype);
        }

//        zygosityCount = new ZygosityCount();
//        for (String sample : getSampleNames()) {
//            Genotype genotype = getGenotype(sample);
//            zygosityCount.incrementCount(genotype);
//        }

        String afString = null;
        String[] alleleFreqKeys = {"AF", "GMAF"};
        try {
            for (String alleleFreqKey : alleleFreqKeys) {
                afString = variantContext.getAttributeAsString(alleleFreqKey, "-1");
                alleleFreqs = parseDoubleArrayString(afString);
                if (alleleFreqs[0] >= 0) break;
            }
        } catch (NumberFormatException e) {
            alleleFreqs = new double[]{-1};
            log.error("Error parsing allele frequency: " + afString);
        }

        String acKey = "AC";
        String acString = variantContext.getAttributeAsString(acKey, null);
        if (acString != null) {
            try {
                alleleCounts = parseIntArrayString(acString);
            } catch (NumberFormatException e) {
                log.error("Error parsing allele counts:" + acString);
            }
        }

        String anKey = "AN";
        String anString = variantContext.getAttributeAsString(anKey, null);
        if(anString != null) {
            try {
                totalAlleleCount = Integer.parseInt(anString);
            } catch (NumberFormatException e) {
                log.error("Error parsing 'AN' attribute: " + anString);
            }
        }

        // The variant is considered "NON_REF" if all of the alternate alleles are non-ref
        for(htsjdk.variant.variantcontext.Allele a : variantContext.getAlternateAlleles()){
            if(a.isNonRefAllele()) {
                nonRef = true;
            } else {
                nonRef = false;
                break;
            }
        }

    }


    /**
     * Allele frequency is a comma separated list of doubles
     * We strip away brackets and parentheses
     *
     * @param afString
     * @return
     */
    private double[] parseDoubleArrayString(String afString) throws NumberFormatException {
        afString = afString.replaceAll("[\\[\\]\\(\\)]", "");
        String[] tokens = afString.split(",");
        double[] result = new double[tokens.length];
        for (int ii = 0; ii < tokens.length; ii++) {
            result[ii] = Double.parseDouble(tokens[ii].trim());
        }
        return result;
    }

    private int[] parseIntArrayString(String afString) throws NumberFormatException {
        afString = afString.replaceAll("[\\[\\]\\(\\)]", "");
        String[] tokens = afString.split(",");
        int[] result = new int[tokens.length];
        for (int ii = 0; ii < tokens.length; ii++) {
            result[ii] = Integer.parseInt(tokens[ii].trim());
        }
        return result;
    }

    /**
     * Compute the average methylation rate for those samples with data (i.e. with methylation rate recorded).
     */
    private void computeMethylationRate() {

        double methTotal = 0;
        int samplesWithData = 0;
        final int size = getSampleNames().size();
        if (size > 0) {
            for (String sample : getSampleNames()) {
                Genotype genotype = getGenotype(sample);
                double mr = genotype.getAttributeAsDouble("MR");
                double goodBaseCount = genotype.getAttributeAsDouble("MR");
                if (!Double.isNaN(mr) && !Double.isNaN(goodBaseCount) && goodBaseCount > VariantTrack.METHYLATION_MIN_BASE_COUNT) {
                    methTotal += mr;
                    samplesWithData++;
                }
            }
            methylationRate = samplesWithData == 0 ? 0 : methTotal / samplesWithData;
            coveredSampleFraction = ((double) samplesWithData) / size;
        }
    }


    public String getID() {
        return variantContext.getID();
    }

    public boolean isFiltered() {
        return variantContext.isFiltered();
    }

    public String getAttributeAsString(String key) {
        return variantContext.getAttributeAsString(key, null);
    }

    public String getReference() {
        return variantContext.getReference().toString();
    }

    public List<Allele> getAlternateAlleles() {
        if (alternateAlleles == null) {
            List<htsjdk.variant.variantcontext.Allele> tmp = variantContext.getAlternateAlleles();
            alternateAlleles = new ArrayList<Allele>(tmp.size());
            for (htsjdk.variant.variantcontext.Allele a : tmp) {
                alternateAlleles.add(new VCFAllele(a));
            }
        }
        return alternateAlleles;
    }

    public double getPhredScaledQual() {
        return variantContext.getPhredScaledQual();
    }

    @Override
    public boolean hasLog10PError() {
        return variantContext.hasLog10PError();
    }

    public String getType() {
        return variantContext.getType().toString();
    }


    public boolean isNonRef() {
        return nonRef;
    }

    /**
     * Return the allele frequency as annotated with an AF or GMAF attribute.  A value of -1 indicates
     * no annotation (unknown allele frequency).
     */
    public double[] getAlleleFreqs() {
        return alleleFreqs;
    }


    @Override
    public double getAlternateAlleleFrequency() {
        double af = 0;
        double[] afreqs = getAlleleFreqs();
        if (afreqs != null) {
            for (int i = 0; i < afreqs.length; i++) {
                af += afreqs[i];
            }
        }
        return af;
    }


    @Override
    public int[] getAlleleCounts() {
        return alleleCounts;
    }

    public int getTotalAlleleCount() {
        return totalAlleleCount;
    }

    public double getAlleleFraction() {
        if(alleleCounts != null && alleleCounts.length > 0 && totalAlleleCount > 0) {
            double ac = 0;
            for(int i=0; i<alleleCounts.length; i++) {
                ac += alleleCounts[i];
            }
            return ac / totalAlleleCount;
        }
        else {
            return -1;
        }
    }


    /**
     * Return the methylation rate as annoted with a MR attribute.  A value of -1 indicates
     * no annotation (unknown methylation rate).  This option is only applicable for dna methylation data.
     */
    public double getMethlationRate() {
        if (Double.isNaN(methylationRate)) {
            computeMethylationRate();
        }
        return methylationRate;
    }

    public double getCoveredSampleFraction() {
        if (Double.isNaN(coveredSampleFraction)) {
            computeMethylationRate();
        }
        return coveredSampleFraction;
    }

    public Collection<String> getSampleNames() {
        return variantContext.getSampleNames();
    }

    public Map<String, Object> getAttributes() {
        return variantContext.getAttributes();
    }

    @Override
    public Genotype getGenotype(String sample) {
        return genotypeMap.get(sample);
    }

    public Collection<String> getFilters() {
        return variantContext.getFilters();
    }

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
        if (this.start < 0) {
            calcStart();
        }
        return this.start;
    }

    @Override
    public int getEnd() {
        String chr2 = variantContext.getAttributeAsString("CHR2", null);
        if(chr2 != null) {
            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            chr2 = genome == null ? chr2 : genome.getCanonicalChrName(chr2);
        }
        if(chr2 == null || chr2.equals(getChr())) {
            return variantContext.getEnd();
        }
        else {
            return this.start + 1;   // An inter-chr variant
        }
    }

    @Override
    public String toString() {
        return String.format("VCFVariant[%s:%d-%d]", getChr(), getStart(), getEnd());
    }

    @Override
    public String getPositionString() {
        if (variantContext.getStart() == variantContext.getEnd()) {
            return numFormat.format(variantContext.getStart());
        } else {
            String chr2 = variantContext.getAttributeAsString("CHR2", null);
            if(chr2 == null || chr2.equals(getChr())) {
             return   numFormat.format( variantContext.getStart()) + "-" + numFormat.format(variantContext.getEnd());
            }
            else {
                // TODO -- I'm not sure this is correct
                return getChr() + ":" + numFormat.format(variantContext.getStart()) + " " +
                        chr2 + ":"  + numFormat.format(variantContext.getEnd());
            }
        }

    }


    public String getSource() {
        return variantContext.getSource();
    }

    public VariantContext getVariantContext() {
        return variantContext;
    }

    /**
     * VCFs specify padding bases at the beginning of indels so they can be positioned properly.
     * We display the variant only where it actually differs from the reference. So we find the longest
     * common prefix between the reference and variants
     */
    private void calcStart() {
        int prefixLength = 0;

        if (variantContext.getType() == VariantContext.Type.INDEL || variantContext.getType() == VariantContext.Type.MIXED) {
            prefixLength = findCommonPrefixLength();
        }

        /**
         * The VCF spec defines POS as the base preceding the polymorphism for non-ref symbolic alleles.
         *
         */
        if (variantContext.getType(true) == VariantContext.Type.SYMBOLIC) {
            prefixLength = 1;
        }

        this.start = (variantContext.getStart() - 1) + prefixLength;
    }

    /**
     * Find the length of the common prefix between the reference and ALL
     * variant alleles
     *
     * @return
     */
    private int findCommonPrefixLength() {
        String ref = variantContext.getReference().getDisplayString();
        int prefixLength = 0;
        boolean foundmisMatch = false;
        for (int refPos = 0; refPos < ref.length(); refPos++) {
            char refChar = ref.charAt(refPos);
            for (Allele var : getAlternateAlleles()) {
                byte[] varBases = var.getBases();
                if (refPos >= varBases.length || varBases[refPos] != refChar) {
                    foundmisMatch = true;
                    break;
                }
            }
            if (foundmisMatch) {
                break;
            } else {
                prefixLength++;
            }
        }
        return prefixLength;
    }
}
