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

package org.broad.igv.cli_plugin;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Ignore;

/**
 * User: jacob
 * Date: 2012-Aug-20
 */
@Ignore
public class AbstractPluginTest extends AbstractHeadlessTest {

    static boolean haveTool;
    protected static String pluginPath;
    static PluginSpecReader reader;
    static String toolPath;
    static PluginSpecReader.Tool tool;

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();

        reader = PluginSpecReader.create(pluginPath);
        initTool(reader);
    }

    /**
     * Sets static variables, and checks that tool exists
     * on this machine
     *
     * @param reader
     */
    public static void initTool(PluginSpecReader reader) {
        tool = reader.getTools().get(0);
        toolPath = reader.getToolPath(tool);
        haveTool = PluginSpecReader.isToolPathValid(toolPath);
        Assume.assumeTrue(haveTool);
    }

    public static PluginSpecReader getCatReader() {
        String path = TestUtils.DATA_DIR + "cli_plugin/cat_plugin.xml";
        PluginSpecReader reader = PluginSpecReader.create(path);
        return reader;
    }

    public static PluginSpecReader.Command findCommandElementByName(PluginSpecReader.Tool tool, String cmd){

        //Find the command element
        PluginSpecReader.Command command = null;
        for (PluginSpecReader.Command curCmd : tool.commandList) {
            if (curCmd.cmd.equalsIgnoreCase(cmd)) {
                return curCmd;
            }
        }
        return null;

    }
}
