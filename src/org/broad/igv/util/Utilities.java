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

package org.broad.igv.util;

//~--- non-JDK imports --------------------------------------------------------

import biz.source_code.base64Coder.Base64Coder;
import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.zip.CRC32;

/**
 *
 */
public class Utilities {

    private static Logger log = Logger.getLogger(Utilities.class);

    public static String base64Encode(String str) {
        return Base64Coder.encodeString(str);

    }

    public static String base64Decode(String str) {
        return Base64Coder.decodeString(str);
    }


    /**
     * Changes a files extension.
     */
    final public static File changeFileExtension(File file, String extension) {

        if ((file == null) || (extension == null) || extension.trim().equals("")) {
            return null;
        }

        String path = file.getAbsolutePath();
        String newPath = "";
        String filename = file.getName().trim();

        if ((filename != null) && !filename.equals("")) {

            int periodIndex = path.lastIndexOf(".");

            newPath = path.substring(0, periodIndex);
            newPath += extension;
        }

        return new File(newPath);
    }

    /**
     * Reads an xml from an input stream and creates DOM document.
     *
     * @param inStream
     * @return
     * @throws ParserConfigurationException
     * @throws IOException
     * @throws SAXException
     */
    public static Document createDOMDocumentFromXmlStream(InputStream inStream)
            throws ParserConfigurationException, IOException, SAXException {

        DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
        Document xmlDocument = documentBuilder.parse(inStream);

        return xmlDocument;
    }

    public static Map<String, String> getAttributes(Node node) {
        HashMap<String, String> attributes = new HashMap();

        NamedNodeMap tNodeMap = node.getAttributes();
        if(tNodeMap != null){
            for (int i = 0; i < tNodeMap.getLength(); i++) {
                Node nd = tNodeMap.item(i);
                String value = nd.getNodeValue();
                if (value != null && value.length() > 0) {
                    attributes.put(nd.getNodeName(), value);
                }
            }
        }
        return attributes;
    }

    public static String getString(Document document) {
        StreamResult streamResult;
        try {

            TransformerFactory factory = TransformerFactory.newInstance();
            Transformer transformer = factory.newTransformer();

            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");

            streamResult = new StreamResult(new StringWriter());
            DOMSource source = new DOMSource(document);
            transformer.transform(source, streamResult);
        } catch (TransformerException e) {
            return null;
        }

        return streamResult.getWriter().toString();
    }


    /**
     * Get the node value of the attribute {@code key} from
     * the set {@code attr}. If that attribute is not present, null is returned
     *
     * @param attr
     * @param key
     * @return
     */
    public static String getNullSafe(NamedNodeMap attr, String key) {
        Node node = attr.getNamedItem(key);
        return node != null ? node.getNodeValue() : null;
    }

    static public Comparator getNumericStringComparator() {

        Comparator comparator = new Comparator<String>() {

            public int compare(String s1, String s2) {
                StringTokenizer st1 = new StringTokenizer(s1, " ");
                StringTokenizer st2 = new StringTokenizer(s2, " ");

                while (st1.hasMoreTokens() && st2.hasMoreTokens()) {
                    String t1 = st1.nextToken();
                    String t2 = st2.nextToken();

                    int c;

                    try {
                        Integer i1 = new Integer(t1);
                        Integer i2 = new Integer(t2);

                        c = i1.compareTo(i2);
                    } catch (NumberFormatException e) {
                        c = t1.compareTo(t2);
                    }

                    if (c != 0) {
                        return c;
                    }
                }

                return 0;
            }
        };

        return comparator;
    }

    /**
     * @param buffer
     * @return
     * @throws java.io.IOException
     */
    static private long getCrc(byte[] buffer) throws IOException {

        CRC32 crc = new CRC32();

        crc.reset();
        crc.update(buffer, 0, buffer.length);
        return crc.getValue();
    }

    /**
     * Parses the file name from a URL.
     *
     * @param url
     * @return
     */
    static public String getFileNameFromURL(String url) {

        int lastIndexOfSlash = url.lastIndexOf("/");

        String fileName = null;

        if (lastIndexOfSlash > -1) {
            fileName = url.substring(lastIndexOfSlash+1, url.length());
        } else {
            fileName = url;
        }

        return fileName;
    }

    public static boolean objectEqual(Object str1, Object str2) {

        if(str1 == null) {
            return str2 == null;
        }
        else if(str2 == null) {
            return str1 == null;
        }
        else {
            return str1.equals(str2);
        }

    }
}
