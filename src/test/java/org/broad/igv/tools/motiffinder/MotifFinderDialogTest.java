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

package org.broad.igv.tools.motiffinder;

import org.broad.igv.AbstractHeadlessTest;
import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * @author jacob
 * @date 2013-Apr-10
 */
public class MotifFinderDialogTest extends AbstractHeadlessTest {

    @Test
    public void testBasicChars() throws Exception{
        String chars = "ACTGN";

        assertTrue(MotifFinderDialog.checkIUPACPatternValid(chars));
        assertTrue(MotifFinderDialog.checkNucleotideRegex(chars));
    }

    @Test
    public void testInvalidChars() throws Exception{
        String chars = "XQZACTGN";

        assertFalse(MotifFinderDialog.checkIUPACPatternValid(chars));
        assertFalse(MotifFinderDialog.checkNucleotideRegex(chars));
    }

    @Test
    public void testAmbiguity() throws Exception{
        String chars = "GATCRYMKSWHBVDN";

        assertTrue(MotifFinderDialog.checkIUPACPatternValid(chars));

        assertFalse(MotifFinderDialog.checkNucleotideRegex(chars));
    }

    @Test
    public void testRegexGood() throws Exception{
        String chars = "TATTAAATTGC[A,T,GCTT]*AA+CGCT{4,5}";

        assertFalse(MotifFinderDialog.checkIUPACPatternValid(chars));

        assertTrue(MotifFinderDialog.checkNucleotideRegex(chars));
    }

    @Test
    public void testRegexBad() throws Exception{
        //R is not a valid letter
        String chars = "TATTAAA{3,10}RN";

        assertFalse(MotifFinderDialog.checkIUPACPatternValid(chars));

        assertFalse(MotifFinderDialog.checkNucleotideRegex(chars));
    }
}
