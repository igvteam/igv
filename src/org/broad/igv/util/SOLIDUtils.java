/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
