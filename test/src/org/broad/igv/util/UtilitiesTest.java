package org.broad.igv.util;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by jrobinso on 12/17/16.
 */
public class UtilitiesTest {
    @Test
    public void isAllLetters() throws Exception {

        boolean b = Utilities.isAllLetters("StringWithNothingButLetters");
        assertTrue(b);

        b = Utilities.isAllLetters("StringWithNumbers123");
        assertFalse(b);

        b = Utilities.isAllLetters("String with spaces");
        assertFalse(b);
    }

}