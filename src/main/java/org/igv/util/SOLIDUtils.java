package org.igv.util;

/**
 * User: jrobinso
 * Date: Feb 21, 2010
 */
public class SOLIDUtils {
    /**
     * Return the sequence in Color Space (SOLID alignment encoding)
     */
    public static byte[] convertToColorSpace(byte[] baseSequence) {
        // We need to know the base just to the left of the start

        if (baseSequence == null || baseSequence.length == 0) {
            return baseSequence;
        }

        byte[] csSequence = new byte[baseSequence.length];
        // Can't know the first base
        csSequence[0] = baseSequence[0]; //
        int c1 = baseToCS(baseSequence[0]);
        for (int i = 1; i < baseSequence.length; i++) {
            int c2 = baseToCS(baseSequence[i]);
            if (c2 == 4 && c1 != 4) {
              csSequence[i] = 4;
            } else if (c1 == 4 && c2 != 4) {
              csSequence[i] = 5;
            } else if (c1 == 4 && c2 == 4) {
             csSequence[i] = 6;
            } else {
                csSequence[i] = (byte) (c1 ^ c2);
            }
            c1 = c2;
        }
        return csSequence;

    }

    public static int baseToCS(byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return 0;
            case 'C':
            case 'c':
                return 1;
            case 'G':
            case 'g':
                return 2;
            case 'T':
            case 't':
                return 3;
            case 'N':
            case 'n':
                return 4;
        }
        return -1;
    }
}
