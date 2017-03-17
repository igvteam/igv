package org.broad.igv.sam.lite;

/**
 * Created by jrobinso on 3/14/17.
 */
class CigarOperator {

    public char letter;
    public int length;

    public CigarOperator(int length, char letter) {
        this.length = length;
        this.letter = letter;
    }
}
