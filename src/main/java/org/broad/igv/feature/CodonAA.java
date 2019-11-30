package org.broad.igv.feature;

public class CodonAA {
    String codon;
    AminoAcid aminoAcid;

    public CodonAA(String codon, AminoAcid aminoAcid) {
        this.codon = codon;
        this.aminoAcid = aminoAcid;
    }

    public String getCodon() {
        return codon;
    }

    public AminoAcid getAminoAcid() {
        return aminoAcid;
    }

    public char getSymbol() {
        return aminoAcid.getSymbol();
    }
}
