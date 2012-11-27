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

package org.broad.igv.dev.plugin;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.w3c.dom.Element;

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
    static Element tool;

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();

        reader = PluginSpecReader.create(pluginPath);
        tool = reader.getTools().get(0);
        toolPath = tool.getAttribute("default_path");
        haveTool = PluginSpecReader.isToolPathValid(toolPath);
        Assume.assumeTrue(haveTool);
    }

    public static PluginSpecReader getCatReader(){
        String path = TestUtils.DATA_DIR + "plugin/cat_plugin.xml";
        PluginSpecReader reader = PluginSpecReader.create(path);
        return reader;
    }
}
