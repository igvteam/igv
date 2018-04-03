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
import org.broad.igv.Globals;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assume;
import org.junit.Test;

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
        PluginFeatureSource pluginSource = new PluginFeatureSource(cmd, arguments, command.outputList.get(0), reader.getSpecPath());

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
