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

package org.broad.igv.ui;

import org.junit.After;
import org.junit.Assume;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012/02/27
 */
public class ReadmeParserTest {

    private String path = "docs/igvtools_readme.txt";
    private ReadmeParser parser;

    @Before
    public void setUp() {
        parser = new ReadmeParser(path);
    }

    @After
    public void tearDown() {
        parser = null;
    }

    @Test
    public void testGetBadCmd() throws Exception {
        String fake_cmd = "fake_cmd";
        String info = parser.getDocForCommand(fake_cmd);
        assertTrue(info.contains("Command " + fake_cmd + " not found"));
    }


    @Test
    public void testGetCommands() throws Exception {
        String[] commands = {"count", "sort", "tile", "index", "sort", "version", "toTDF", "formatexp", "gui"};
        for (String cmd : commands) {
            tstGetCommand(cmd);
        }

    }

    public void tstGetCommand(String cmd) throws Exception {
        String info = parser.getDocForCommand(cmd);
        assertNotNull(info);
        assertFalse(info.contains("Command " + cmd + " not found"));
        assertFalse(info.matches("Command .* not found.*"));
        assertTrue(info.length() > 0);
    }

    //This test is here to test performance of Jenkins on test failure
    @Test
    public void setTestFailByCommandLine() throws Exception {
        String makefail = System.getProperty("make.fail", "false");
        boolean shouldFail = Boolean.parseBoolean(makefail);
        //Test ignored unless we specify failure, because it's not really a test
        Assume.assumeTrue(shouldFail);

        assertFalse(shouldFail);
    }
}
