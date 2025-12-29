/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.aa;

/**
 * @author jrobinso
 */
public class AminoAcid {

    private String fullName;
    private String abbrevName;
    private char symbol;

    public static AminoAcid NULL_AMINO_ACID = new AminoAcid("", "", ' ');

    AminoAcid(String fullName, String abbrevName, char symbol) {
        this.fullName = fullName;
        this.abbrevName = abbrevName;
        this.symbol = symbol;
    }

    public String getFullName() {
        return fullName;
    }

    public String getShortName() {
        return abbrevName;
    }

    public char getSymbol() {
        return symbol;
    }

    public boolean equalsByName(String mutAA) {
        return fullName.equals(mutAA) || abbrevName.equals(mutAA)
                || String.valueOf(symbol).equals(mutAA);
    }
}
