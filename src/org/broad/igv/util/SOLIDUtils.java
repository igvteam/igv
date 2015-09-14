/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.util;

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
