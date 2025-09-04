package org.broad.igv.sam.mods;

import htsjdk.samtools.util.SequenceUtil;

import java.util.HashMap;
import java.util.Map;

/**
 * A set of base modification likelihoods for a given base in a read
 *
 * @author jrobinso
 * @date Jul 2024
 */
public class BaseModificationSet {

    char base;
    char strand;
    String modification;
    Map<Integer, Byte> likelihoods;
    char canonicalBase;

    public BaseModificationSet(char base, char strand, String modification,  Map<Integer, Byte> likelihoods) {
        this.base = base;
        this.modification = modification;
        this.strand = strand;
        this.likelihoods = likelihoods;
        this.canonicalBase = strand == '+' ? base : (char) SequenceUtil.complement((byte) base);
    }

    public char getBase() {
        return base;
    }

    public char getCanonicalBase() {
        return canonicalBase;
    }

    public String getModification() {
        return modification;
    }

    public char getStrand() {
        return strand;
    }

    public Map<Integer, Byte> getLikelihoods() {
        return likelihoods;
    }

    public boolean containsPosition(Integer pos) {
        return likelihoods.containsKey(pos);
    }

    /**
     * Return a descriptive string for the modification at the given position of the read sequence*
     * @param pos - position in the read sequence  (left to right, as recorded in BAM record, not 5'->3')
     * @return
     */
    public String valueString(int pos) {
        int l = (int) (100.0 * Byte.toUnsignedInt(likelihoods.get(pos)) / 255);
        return "Base modification: " +
                ((codeValues.containsKey(modification)) ? codeValues.get(modification) : modification) +  " (" + l + "%)";
    }

    static Map<String, String> codeValues;

    static {
        codeValues = new HashMap<>();
        codeValues.put("m", "5mC");
        codeValues.put("h", "5hmC");
        codeValues.put("f", "5fC");
        codeValues.put("c", "5caC");
        codeValues.put("g", "5hmU");
        codeValues.put("e", "5fU");
        codeValues.put("b", "5caU");
        codeValues.put("a", "6mA");
        codeValues.put("o", "8xoG");
        codeValues.put("n", "Xao");
        codeValues.put("C", "Unknown C");
        codeValues.put("T", "Unknown T");
        codeValues.put("A", "Unknown A");
        codeValues.put("G", "Unknown G");
        codeValues.put("N", "Unknown");

    }
}
