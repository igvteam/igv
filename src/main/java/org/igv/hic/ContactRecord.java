package org.igv.hic;

public record ContactRecord(int bin1, int bin2, float counts, float normCounts) {

    public ContactRecord(int bin1, int bin2, float counts) {
        this(bin1, bin2, counts, counts);
    }

    public String getKey() {
        return bin1 + "_" + bin2;
    }
}