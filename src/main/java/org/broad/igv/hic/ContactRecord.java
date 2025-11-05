package org.broad.igv.hic;

public record ContactRecord(int bin1, int bin2, double counts) {
    public String getKey() {
        return bin1 + "_" + bin2;
    }
}