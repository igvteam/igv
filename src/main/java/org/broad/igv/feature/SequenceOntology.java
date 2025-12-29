package org.broad.igv.feature;

import java.util.HashSet;
import java.util.Set;

/**
 * A static class for defining identifiers from the Sequence Ontology project and related efforts around GFF, EMBL, and
 * NCBI formats.
 * <p/>
 * Reference:    http://www.sequenceontology.org/
 */
public class SequenceOntology {


    public static Set<String> fivePrimeUTRTypes = new HashSet<String>();
    public static Set<String> threePrimeUTRTypes = new HashSet<String>();
    public static Set<String> utrTypes = new HashSet<String>();
    public static Set<String> cdsTypes = new HashSet<String>();
    public static Set<String> exonTypes = new HashSet<String>();
    public static Set<String> transcriptParts = new HashSet();
    public static Set<String> mrnaParts = new HashSet();
    public static Set<String> geneParts = new HashSet();


    static {
        fivePrimeUTRTypes.add("five_prime_UTR");
        fivePrimeUTRTypes.add("5'-UTR");
        fivePrimeUTRTypes.add("5'-utr");
        fivePrimeUTRTypes.add("5UTR");

        threePrimeUTRTypes.add("three_prime_UTR");
        threePrimeUTRTypes.add("3'-utr");
        threePrimeUTRTypes.add("3'-UTR");
        threePrimeUTRTypes.add("3UTR");

        utrTypes.addAll(SequenceOntology.fivePrimeUTRTypes);
        utrTypes.addAll(SequenceOntology.threePrimeUTRTypes);
        utrTypes.add("UTR");    // Some Ensemble gtf files
        cdsTypes.add("CDS");
        cdsTypes.add("cds");

        exonTypes.add("exon");
        exonTypes.add("coding_exon");

        mrnaParts.addAll(utrTypes);
        mrnaParts.addAll(cdsTypes);
        mrnaParts.addAll(exonTypes);
        mrnaParts.add("intron");
        mrnaParts.add("polyA_sequence");
        mrnaParts.add("polyA_site");
        mrnaParts.add("start_codon");
        mrnaParts.add("stop_codon");

        transcriptParts.addAll(mrnaParts);

        geneParts.addAll(transcriptParts);
        geneParts.add("transcript");
        geneParts.add("processed_transcript");
        geneParts.add("mrna");
        geneParts.add("mRNA");

    }

    public static boolean isCoding(String type) {
        return cdsTypes.contains(type) || type.equals("coding_exon");
    }
}
