package org.igv.hic;

public record ContactRecord(int bin1, int bin2, float counts) {
    public String getKey() {
        return bin1 + "_" + bin2;
    }
}