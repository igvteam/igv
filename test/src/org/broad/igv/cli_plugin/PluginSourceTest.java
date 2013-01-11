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

package org.broad.igv.cli_plugin;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assume;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.namespace.QName;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * User: jacob
 * Date: 2013-Jan-03
 */
public class PluginSourceTest extends AbstractHeadlessTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 30e6);

    @Test
    public void testMarshall() throws Exception{
        Assume.assumeTrue(!Globals.IS_WINDOWS);

        PluginSpecReader reader = AbstractPluginTest.getCatReader();

        PluginSpecReader.Tool tool = reader.getTools().get(0);
        PluginSpecReader.Command command = tool.commandList.get(0);
        List<Argument> argumentList = command.argumentList;


        LinkedHashMap<Argument, Object> arguments = new LinkedHashMap<Argument, Object>(argumentList.size());

        int argnum = 0;
        arguments.put(argumentList.get(argnum++), "");

        TrackLoader loader = new TrackLoader();
        String[] paths = new String[]{TestUtils.DATA_DIR + "bed/test.bed", TestUtils.DATA_DIR + "bed/testAlternateColor.bed"};
        for (String path : paths) {
            TestUtils.createIndex(path);
            FeatureTrack track = (FeatureTrack) loader.load(new ResourceLocator(path), genome).get(0);
            arguments.put(argumentList.get(argnum++), track);
        }


        List<String> cmd = Arrays.asList(reader.getToolPath(tool), command.cmd);
        PluginFeatureSource pluginSource = new PluginFeatureSource(cmd, arguments, command.parser, reader.getSpecPath());

        //-------------------------//

//        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
//        DocumentBuilder builder = factory.newDocumentBuilder();
//        Document doc = builder.newDocument();

        JAXBContext jc = JAXBContext.newInstance(PluginSource.class);
        Marshaller m = jc.createMarshaller();
        m.setProperty(Marshaller.JAXB_FRAGMENT, true);

        JAXBElement inel = new JAXBElement(new QName("", "source"), PluginSource.class, pluginSource);
        m.marshal(inel, System.out);
    }

//    @Test
//    public void testHashMapMarshall() throws Exception{
//        JAXBContext jc = JAXBContext.newInstance(MyTestClass.class);
//        Marshaller m = jc.createMarshaller();
//        m.setProperty(Marshaller.JAXB_FRAGMENT, true);
//
//        MyTestClass testObj = new MyTestClass();
//        testObj.myMap.put("ka", "va");
//        testObj.myMap.put("kb", "vb");
//
//        JAXBElement inel = new JAXBElement(new QName("", "source"), MyTestClass.class, testObj);
//        m.marshal(inel, System.out);
//    }
//
//    static class MyTestClass{
//
//        @XmlJavaTypeAdapter(MyMapAdapter.class)
//        public LinkedHashMap<String, String> myMap = new LinkedHashMap<String, String>();
//
//    }
//
//    static class XmlMap{
//        public List<XmlMapEntry> entry =
//                new ArrayList<XmlMapEntry>();
//    }
//
//    static class XmlMapEntry{
//        @XmlAttribute
//        public String key;
//
//        @XmlValue
//        public String value;
//    }
//
//    public static final class MyMapAdapter extends XmlAdapter<XmlMap, Map<String, String>> {
//
//        @Override
//        public LinkedHashMap<String, String> unmarshal(XmlMap v) throws Exception {
//            return null; //TODO
//        }
//
//        @Override
//        public XmlMap marshal(Map<String, String> v) throws Exception {
//            XmlMap xmlMap = new XmlMap();
//            for(Map.Entry<String, String> loopEntry: v.entrySet()){
//                XmlMapEntry entry = new XmlMapEntry();
//                entry.key = loopEntry.getKey();
//                entry.value = loopEntry.getValue();
//                xmlMap.entry.add(entry);
//            }
//            return xmlMap;
//        }
//    }

}
