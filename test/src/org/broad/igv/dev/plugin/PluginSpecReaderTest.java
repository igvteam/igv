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

import org.broad.igv.PreferenceManager;
import org.broad.igv.util.TestUtils;
import org.junit.Test;
import org.w3c.dom.Element;

import java.io.File;
import java.io.FilenameFilter;
import java.lang.reflect.Field;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012-Aug-09
 */
public class PluginSpecReaderTest {

    @Test
    public void testLoadExternalClassFromJar() throws Exception {
        File extDir = new File(TestUtils.DATA_DIR + "plugin", "StubJar.jar");
        tstLoadExternalClass(extDir);
    }

    @Test
    public void testLoadExternalClassFromDir() throws Exception {
        File extDir = new File(TestUtils.DATA_DIR + "plugin", "/");
        tstLoadExternalClass(extDir);
    }

    public void tstLoadExternalClass(File extDir) throws Exception {

        URL homeURL = new URL("file:" + extDir.getAbsolutePath());

        ClassLoader loader = URLClassLoader.newInstance(
                new URL[]{homeURL},
                getClass().getClassLoader()
        );

        Class clazz = loader.loadClass("org.broad.igv.StubClass");
        Object obj = clazz.newInstance();

        Field name = clazz.getDeclaredField("name");

        String val = (String) name.get(obj);
        assertTrue(val.contains("I am a stub class"));
    }

    @Test
    public void testReadSpec() throws Exception {
        String path = TestUtils.DATA_DIR + "plugin/cat_plugin.xml";
        PluginSpecReader reader = PluginSpecReader.create(path);
        assertNotNull(reader.document);

        List<Element> tools = reader.getTools();
        assertEquals(1, tools.size());

        List<Element> commands = reader.getCommands(tools.get(0));
        assertEquals(1, commands.size());

        List<Argument> arguments = reader.getArguments(tools.get(0), commands.get(0));
        assertEquals(3, arguments.size());
    }

    /**
     * Check that we fail fast when reading an XML document which is NOT a plugin reader
     *
     * @throws Exception
     */
    @Test
    public void testReadSpecFail() throws Exception {
        String path = TestUtils.DATA_DIR + "sessions/testBedsRelPath.xml";
        PluginSpecReader reader = PluginSpecReader.create(path);
        assertNull(reader);
    }

    @Test
    public void testCustomToolPath() throws Exception{
        String path = TestUtils.DATA_DIR + "plugin/cat_plugin.xml";
        PluginSpecReader reader = PluginSpecReader.create(path);
        Element tool = reader.getTools().get(0);
        String toolName = tool.getAttribute(PluginSpecReader.TOOL_NAME_KEY);
        assertEquals("cat", toolName);

        String newpath = "/dev/zero";
        PreferenceManager.getInstance().putPluginPath(reader.getId(), toolName, newpath);
        assertEquals(newpath, reader.getToolPath(tool));
    }


    //Check that all of the plugins we specify exist
    @Test
    public void testBuiltinPluginsValid() throws Exception {
        List<String> pluginNames = PluginSpecReader.getBuiltinPlugins();
        for (String pluginName : pluginNames) {
            String relPath = "resources/" + pluginName;
            URL url = PluginSpecReader.class.getResource(relPath);
            //System.out.println(url);
            assertNotNull(url);
        }
    }

    //Check that each plugin file is in the contents file
    @Test
    public void testBuiltinPluginsComplete() throws Exception {
        String pluginsPath = "src/" + PluginSpecReader.class.getPackage().getName().replace('.', '/');
        File pluginResourceDir = new File(pluginsPath, "resources");
        String[] fileNames = pluginResourceDir.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".xml");
            }
        });

        Set<String> actPlugins = new HashSet<String>(Arrays.asList(fileNames));
        Set<String> expPlugins = new HashSet<String>(PluginSpecReader.getBuiltinPlugins());

        assertEquals(expPlugins, actPlugins);

    }


}
