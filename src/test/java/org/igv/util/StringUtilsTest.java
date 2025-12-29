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
