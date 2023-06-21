package org.broad.igv.sam.mods;

import htsjdk.samtools.util.SequenceUtil;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

public class BaseModificationKey {
    char base;
    char strand;
    String modification;

    static Map<String, BaseModificationKey> keyCache = Collections.synchronizedMap(new HashMap<>());

    public static BaseModificationKey getKey(char base, char strand, String modification) {

        String name = base + strand + modification;
        if (keyCache.containsKey(name)) {
            return keyCache.get(name);
        } else {
            BaseModificationKey key = new BaseModificationKey(base, strand, modification);
            keyCache.put(name, key);
            return key;
        }

    }

    private BaseModificationKey(char base, char strand, String modification) {
        this.base = base;
        this.strand = strand;
        this.modification = modification;
    }

    public char getBase() {
        return base;
    }

    public char getCanonicalBase() {
        return strand == '+' ? base : (char) SequenceUtil.complement((byte) base);
    }

    public char getStrand() {
        return strand;
    }

    public String getModification() {
        return modification;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        BaseModificationKey key = (BaseModificationKey) o;
        return base == key.base && strand == key.strand && modification.equals(key.modification);
    }

    @Override
    public int hashCode() {
        return Objects.hash(base, strand, modification);
    }

    @Override
    public String toString() {
        return "" + base + strand + modification;
    }
}
