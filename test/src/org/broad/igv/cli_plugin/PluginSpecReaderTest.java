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

import org.broad.igv.PreferenceManager;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

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
        File extDir = new File(TestUtils.DATA_DIR + "cli_plugin", "StubJar.jar");
        tstLoadExternalClass(extDir);
    }

    @Test
    public void testLoadExternalClassFromDir() throws Exception {
        File extDir = new File(TestUtils.DATA_DIR + "cli_plugin", "/");
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
    public void testReadCatSpec() throws Exception {
        PluginSpecReader reader = AbstractPluginTest.getCatReader();
        assertNotNull(reader.document);

        List<PluginSpecReader.Tool> tools = reader.getTools();
        assertEquals(1, tools.size());

        List<PluginSpecReader.Command> commands = tools.get(0).commandList;
        assertEquals(1, commands.size());

        List<Argument> arguments = commands.get(0).argumentList;
        assertEquals(3, arguments.size());

        boolean defOutput = true;
        String[] encCodecs = {null, null, "org.broad.igv.feature.tribble.IGVBEDCodec"};
        for(int ii=0; ii < arguments.size(); ii++){
            Argument arg = arguments.get(ii);
            assertEquals(defOutput, arg.isOutput());
            assertEquals(encCodecs[ii], arg.getEncodingCodec());
        }
    }

    /**
     * Check that we can read parsing attributes
     * @throws Exception
     */
    @Test
    public void testReadParser() throws Exception{
        String path = "resources/bedtools_plugin.xml";
        PluginSpecReader reader = PluginSpecReader.create(path);
        PluginSpecReader.Tool tool = reader.getTools().get(0);
        List<PluginSpecReader.Command> commands = tool.commandList;
        int ind = 0;
        PluginSpecReader.Command multiinter_command = commands.get(ind);
        while(!multiinter_command.cmd.equals("multiinter")){
            multiinter_command = commands.get(ind++);
        }

        PluginSpecReader.Parser parser = multiinter_command.outputList.get(0).parser;
        assertEquals("bed", parser.format);
        assertEquals(true, parser.strict);
        assertTrue(parser.decodingCodec.contains("BEDToolsDecoder"));
    }

    /**
     * Check that we fail fast when reading an XML document which is NOT a cli_plugin reader
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
        String path = TestUtils.DATA_DIR + "cli_plugin/cat_plugin.xml";
        PluginSpecReader reader = PluginSpecReader.create(path);
        PluginSpecReader.Tool tool = reader.getTools().get(0);
        String toolName = tool.name;
        assertEquals("cat", toolName);

        String newpath = "/dev/zero";
        PreferenceManager.getInstance().putToolPath(reader.getId(), toolName, newpath);
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

    //Check that each cli_plugin file is in the contents file
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
