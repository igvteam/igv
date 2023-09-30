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

package org.broad.igv.feature.genome;

/**
 * Created with IntelliJ IDEA.
 * User: jrobinso
 * Date: 6/25/12
 * Time: 11:00 AM
 * To change this template use File | Settings | File Templates.
 */
public class ColorSpaceSequence   {

    Sequence sequence;

    public ColorSpaceSequence(Sequence sequence) {
        this.sequence = sequence;
    }

    /**
     * Return the sequence in Color Space (SOLID alignment encoding)
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public byte[] getSequence(String chr, int start, int end) {
        // We need to know the base just to the left of the start
        int csStart = (start == 0 ? 0 : start - 1);
        byte[] baseSequence = sequence.getSequence(chr, csStart, end);
        if (baseSequence == null || baseSequence.length == 0) {
            return baseSequence;
        }

        byte[] csSequence = new byte[end - start];
        int i = 0;
        int c1 = start == 0 ? 0 : baseToCS(baseSequence[i++]);
        for (; i < baseSequence.length; i++) {
            int c2 = baseToCS(baseSequence[i]);
            csSequence[i] = (byte) (c1 ^ c2);
        }
        return csSequence;

    }

    private static int baseToCS(byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return 0;
            case 'C':
            case 'c':
                return 1;
            case 'T':
            case 't':
                return 2;
            case 'G':
            case 'g':
                return 3;
        }
        return -1;
    }


    public byte getBase(String chr, int position) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
