/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
