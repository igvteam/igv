/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.util;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.StringTokenizer;
import java.util.zip.CRC32;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 *
 */
public class Utilities {

    private static Logger log = Logger.getLogger(Utilities.class);
    final static int ZIP_ENTRY_CHUNK_SIZE = 64000;


    // TODO -- replace use of sun package
    public static String base64Encode(String str) {
        sun.misc.BASE64Encoder encoder = new sun.misc.BASE64Encoder();
        byte[] bytes = str.getBytes();
        return encoder.encode(bytes);

    }

    // TODO -- replace use of sun package
    public static String base64Decode(String str) {
        try {
            return new String((new sun.misc.BASE64Decoder()).decodeBuffer(str));
        } catch (IOException e) {
            log.error("Error decoding string: " + str, e);
            return str;
        }
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
     * Reads an xml from an input file and creates DOM document.
     *
     * @param file
     * @return
     * @throws ParserConfigurationException
     * @throws IOException
     * @throws SAXException
     */
    public static Document createDOMDocumentFromXmlFile(File file)
            throws ParserConfigurationException, IOException, SAXException {

        if (!file.exists()) {
            throw new RuntimeException("Could not find XML file: " + file.getAbsolutePath());
        }

        DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
        Document xmlDocument = documentBuilder.parse(file);

        return xmlDocument;
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
     * Build a ZipEntry.
     *
     * @param entryName
     * @param bytes
     * @param zipOutputStream
     * @param isCompressed
     * @return
     * @throws java.io.IOException
     * @throws java.io.IOException
     */
    public static ZipEntry createZipEntry(String entryName, byte[] bytes,
                                          ZipOutputStream zipOutputStream, boolean isCompressed)
            throws IOException, IOException {

        ZipEntry zipEntry = new ZipEntry(entryName);

        if (isCompressed) {
            zipEntry.setMethod(ZipEntry.DEFLATED);
        } else {
            zipEntry.setMethod(ZipEntry.STORED);
            zipEntry.setCompressedSize(bytes.length);
            zipEntry.setSize(bytes.length);
            zipEntry.setCrc(getCrc(bytes));
        }

        zipOutputStream.putNextEntry(zipEntry);

        BufferedOutputStream bufferedOutputStream = new BufferedOutputStream(zipOutputStream);

        bufferedOutputStream.write(bytes);
        bufferedOutputStream.flush();

        zipOutputStream.closeEntry();

        return zipEntry;
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
            fileName = url.substring(lastIndexOfSlash, url.length());
        } else {
            fileName = url;
        }

        return fileName;
    }

    static public void createZipFile(File zipOutputFile, File[] inputFiles)
            throws FileNotFoundException, IOException {

        if (zipOutputFile == null) {
            return;
        }

        if ((inputFiles == null) || (inputFiles.length == 0)) {
            return;
        }

        ZipOutputStream zipOutputStream = null;

        try {
            zipOutputStream = new ZipOutputStream(new FileOutputStream(zipOutputFile));

            for (File file : inputFiles) {

                if (file == null) {
                    continue;
                }

                long fileLength = file.length();

                ZipEntry zipEntry = new ZipEntry(file.getName());

                zipEntry.setSize(fileLength);
                zipOutputStream.putNextEntry(zipEntry);

                BufferedInputStream bufferedInputstream = null;

                try {
                    InputStream inputStream = new FileInputStream(file);

                    bufferedInputstream = new BufferedInputStream(inputStream);

                    int bytesRead = 0;
                    byte[] data = new byte[ZIP_ENTRY_CHUNK_SIZE];

                    while ((bytesRead = bufferedInputstream.read(data)) != -1) {
                        zipOutputStream.write(data, 0, bytesRead);
                    }
                } finally {
                    if (bufferedInputstream != null) {
                        bufferedInputstream.close();
                    }
                }
            }
        } finally {
            if (zipOutputStream != null) {
                zipOutputStream.flush();
                zipOutputStream.close();
            }
        }
    }

    public static File getFileFromURL(URL url) {
        return getFileFromURLString(url.toExternalForm());
    }

    public static File getFileFromURLString(String urlText) {

        File file = null;

        try {
            URL url = new URL(urlText);
            URI uri = url.toURI();

            file = new File(uri.getPath());
        } catch (MalformedURLException e) {
            log.error(e.getMessage(), e);
        } catch (URISyntaxException e) {
            log.error(e.getMessage(), e);
        }

        return file;
    }

    final public static String[] fastSplit(String text, String delimeter) {

        List<String> list = new ArrayList();

        if (text == null || text.length() < 1) {
            return new String[0];
        }

        while (true) {

            int index = text.indexOf(delimeter);

            // No more text to split
            if (index == -1) {
                if (text.length() > 0) {
                    list.add(text);
                }
                break;
            }

            if (text.length() > 0) {

                // Store text up to index
                String newText = text.substring(0, index);

                // Trim leading delimeters
                while (newText.startsWith(delimeter)) {
                    newText = newText.substring(1, newText.length());
                }
                if (newText.length() > 0) {
                    list.add(newText);
                }

                // Reset the string to the remaining text
                int length = text.length();
                if (length > 1) {
                    text = text.substring(index + 1, length);

                    // Trim leading delimeters
                    while (text.startsWith(delimeter)) {
                        text = text.substring(1, text.length());
                    }
                }
            }
        }

        return list.toArray(new String[list.size()]);
    }
}
