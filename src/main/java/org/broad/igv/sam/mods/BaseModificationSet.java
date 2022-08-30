package org.broad.igv.sam.mods;

import htsjdk.samtools.util.SequenceUtil;

import java.util.Map;

public class BaseModificationSet {

    char base;
    char strand;

    String modification;
    Map<Integer, Byte> likelihoods;

    public BaseModificationSet(char base, char strand, String modification,  Map<Integer, Byte> likelihoods) {
        this.base = base;
        this.modification = modification;
        this.strand = strand;
        this.likelihoods = likelihoods;
    }

    public char getBase() {
        return base;
    }

    public char getCanonicalBase() {
        return strand == '+' ? base : (char) SequenceUtil.complement((byte) base);
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

    public boolean is5mC() {
        return modification.equals("m") &&  ((base == 'C' && strand == '+') || (base == 'G' && strand == '-'));
    }
}
