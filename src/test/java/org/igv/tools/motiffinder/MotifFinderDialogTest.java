package org.igv.tools.motiffinder;

import org.igv.AbstractHeadlessTest;
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
