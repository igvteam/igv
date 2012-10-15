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

import org.broad.igv.util.TestUtils;
import org.junit.Test;
import org.w3c.dom.Element;

import java.io.File;
import java.lang.reflect.Field;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.List;

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


}
