package org.igv.feature.genome;

public class SeqUtils {
    /**
     * Return the complement of a nucleotide or ambiguity code
     *
     * @param inputChar
     * @return
     */
    public static char complementChar(char inputChar) {
        switch (inputChar) {
            case 'A':
                return 'T';
            case 'T':
            case 'U':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
            case 'M':
                return 'K';
            case 'R':
                return 'Y';
            case 'Y':
                return 'R';
            case 'K':
                return 'M';
            case 'V':
                return 'B';
            case 'H':
                return 'D';
            case 'D':
                return 'H';
            case 'B':
                return 'V';

            case 'a':
                return 't';
            case 't':
            case 'u':
                return 'a';
            case 'g':
                return 'c';
            case 'c':
                return 'g';
            case 'm':
                return 'k';
            case 'r':
                return 'y';
            case 'y':
                return 'r';
            case 'k':
                return 'm';
            case 'v':
                return 'b';
            case 'h':
                return 'd';
            case 'd':
                return 'h';
            case 'b':
                return 'v';

            default:
                return inputChar;
        }
    }
}
