package org.igv.util;

import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jun 3, 2010
 * Time: 5:21:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class StringUtilsTest {

    @Test
    public void testBreakQuotedString() throws Exception {

        String quotedString = "abc,def,\"beforeComma,afterComma\"";
        List<String> tokens = StringUtils.breakQuotedString(quotedString, ',');
        assertEquals(3, tokens.size());
        //for(String s : tokens)  System.out.println(s);
        assertEquals(tokens.get(0), "abc");
        assertEquals(tokens.get(1), "def");
        assertEquals(tokens.get(2), "\"beforeComma,afterComma\"");

        quotedString = "load \"test/data/folder with spaces/test,with,commas.wig\"";
        tokens = StringUtils.breakQuotedString(quotedString, ' ');
        assertEquals(2, tokens.size());
        assertEquals("\"test/data/folder with spaces/test,with,commas.wig\"", tokens.get(1));
    }

    @Test
    public void testBreakQuotedStrin2() throws Exception {

        String quotedString = "abc,def,'beforeComma,afterComma'";
        List<String> tokens = StringUtils.breakQuotedString(quotedString, ',');
        assertEquals(3, tokens.size());
        assertEquals(tokens.get(0), "abc");
        assertEquals(tokens.get(1), "def");
        assertEquals(tokens.get(2), "'beforeComma,afterComma'");
    }

    @Test
    public void testIsSmallPositiveInteger() {
        assertTrue(StringUtils.isSmallPositiveInteger("99"));
        assertFalse(StringUtils.isSmallPositiveInteger("100"));
        assertFalse(StringUtils.isSmallPositiveInteger("1.5"));
        assertFalse(StringUtils.isSmallPositiveInteger("-1"));
        assertFalse(StringUtils.isSmallPositiveInteger("foo"));
    }

}
