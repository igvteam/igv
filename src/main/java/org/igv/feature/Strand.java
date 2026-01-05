
package org.igv.feature;


public enum Strand {

    POSITIVE("+"),
    NEGATIVE("-"),
    NONE(""); //The order here matters because it sets how strand sorts.

    private final String shortString;

    private Strand(String shortString) {
        this.shortString = shortString;
    }

    public static Strand fromString(String strandString) {
        return strandString.equals("+") || strandString.equalsIgnoreCase("POSITIVE")
                ? POSITIVE : (strandString.equals("-") || strandString.equalsIgnoreCase("NEGATIVE")
                ? NEGATIVE : NONE);
    }

    public String toShortString() {
        return shortString;
    }

    public String toString() {
        return this == POSITIVE ? "+" : this == NEGATIVE ? "-" : "";
    }


}
