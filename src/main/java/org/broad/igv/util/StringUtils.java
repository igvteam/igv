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

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.WeakHashMap;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Oct 30, 2009
 * Time: 1:51:36 AM
 * To change this template use File | Settings | File Templates.
 */
public class StringUtils {

    private static Map<String, String> internedStrings = new WeakHashMap<String, String>();

    /**
     * Creates or retrieves an interned copy of {@code string}. This way,
     * we only keep one reference to strings of the same value.
     * Backed by a WeakHashMap
     *
     * @param string
     * @return
     */
    public static String intern(String string) {
        if (!internedStrings.containsKey(string)) {
            internedStrings.put(string, string);
        }
        return internedStrings.get(string);
    }


    public static List<String> breakQuotedString(String string, char splitToken) {
        if(string == null) return null;
        ArrayList<String> strings = new ArrayList<String >();
        if (string.length() == 0) {
            return strings;
        }

        char[] characters = string.toCharArray();
        char c;
        boolean isQuoted = false;
        StringBuffer buff = new StringBuffer(100);
        for (int i = 0; i < characters.length; i++) {
            c = characters[i];
            if (isQuoted) {
                if (c == '"') {
                    isQuoted = false;
                }
                buff.append(c);
            } else if (c == '"') {
                isQuoted = true;
                buff.append(c);
            } else {
                if (c == splitToken) {
                    strings.add(buff.toString().trim());
                    buff.setLength(0);
                } else {
                    buff.append(c);
                }
            }
        }
        if (buff.length() > 0) {
            strings.add(buff.toString().trim());
        }
        return strings;

    }

    /**
     * Return a possibly shortened representation of the input string
     */
    public static String checkLength(String string, int maxLength) {

        if (string.length() <= maxLength) {
            return string;
        }

        int nDots = maxLength > 10 ? 3 : (maxLength > 5 ? 2 : 1);
        int m = (Math.max(1, (maxLength - nDots) / 2));
        StringBuffer newString = new StringBuffer(maxLength);
        newString.append(string.substring(0, m));
        for (int i = 0; i < nDots; i++) {
            newString.append('.');
        }
        newString.append(string.substring(string.length() - m));
        return newString.toString();

    }


    /**
     * Covert a genotype string to a 8 byte integer representation
     * <p/>
     * example
     * String genotype = "AC";
     * short genoShort = genoToShort(genotype);
     * <p/>
     * char allel1 = (char) ((genoShort >>> 8) & 0xFF);
     * char allel2 = (char) ((genoShort >>> 0) & 0xFF);
     * <p/>
     * System.out.println("Allele1: " + allel1);
     * System.out.println("Allele2: " + allel2);
     *
     * @param genotype
     * @return
     */
    public static short genoToShort(String genotype) {
        byte[] bytes = genotype.getBytes();
        return (short) ((bytes[0] & 0xff) << 8 | (bytes[1] & 0xff));
    }

    public static String readString(ByteBuffer byteBuffer) throws IOException {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        byte b = -1;
        while ((b = byteBuffer.get()) != 0) {
            bytes.write(b);
        }
        return new String(bytes.toByteArray());
    }

    /**
     * Decode according to UTF-8. In the extremely unlikely
     * event that we are running on a platform which does not
     * support UTF-8 (it's part of the Java spec), a runtime
     * exception is thrown
     *
     * @param s
     * @return
     */
    public static String decodeURL(String s) {
        if (s == null) {
            return null;
        }
        try {
            return URLDecoder.decode(s, "UTF-8");
        } catch (UnsupportedEncodingException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Encode according to UTF-8. We encode
     * spaces as "%20" rather than "+", seems
     * to work better with others.
     * @param s
     * @return
     */
    public static String encodeURL(String s) {
        if (s == null) {
            return null;
        }
        try {
            s = URLEncoder.encode(s, "UTF-8");
            s = s.replace("+", "%20");
            return s;
        } catch (UnsupportedEncodingException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Attempts to encode only the query portion of a URL,
     * ignoring delimiter characters. If it looks like it's already
     * URL-encoded, does nothing
     * @param url
     * @return
     */
    public static URL encodeURLQueryString(URL url) throws MalformedURLException{

        String[] parts = url.toExternalForm().split("\\?", 2);
        //If there are percent signs we assume it's already url encoded
        if (parts.length > 1 && !parts[1].contains("%")) {
            String queryPart = parts[1];
            String[] params = queryPart.split("&");
            String[] encParms = new String[params.length];
            int ii = 0;
            for (String kv : params) {
                String[] kvs = kv.split("\\=", 2);
                String encString = StringUtils.encodeURL(kvs[0]);
                if(kvs.length == 2){
                    encString = String.format("%s=%s", encString, StringUtils.encodeURL(kvs[1]));
                }
                encParms[ii++] = encString;
            }
            String encQuery = StringUtils.join(encParms, "&");
            String newPath = parts[0] + "?" + encQuery;
            url = new URL(newPath);
        }

        return url;
    }

    /**
     * Copy a text string to the clipboard
     *
     * @param text
     */
    public static void copyTextToClipboard(String text) {
        StringSelection stringSelection = new StringSelection(text);
        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        clipboard.setContents(stringSelection, null);
    }

    /**
     * Converts an input string, of any case, into a series
     * of capitalized words.
     * toCapWords("BOB") -> "Bob"
     * toCapWords("bOb is MY FRiend") -> "Bob Is My Friend"
     *
     * @param text
     * @return
     */
    public static String capWords(String text) {
        String res = "";
        boolean capNext = true;
        String s;

        for (char c : text.toLowerCase().toCharArray()) {
            s = Character.toString(c);
            if (capNext) {
                s = s.toUpperCase();
            }
            res += s;
            capNext = " ".equals(s);
        }
        return res;
    }

    public static String join(List list, String separator) {
        return join(list.toArray(), separator);
    }

    /**
     * Join the objects in {@code array} with the given {@code separator}
     * Useful for making a comma separated list, or URL query string
     * @param array
     * @param separator
     * @return
     */
    public static String join(Object[] array, String separator) {
        int num = array.length;
        if(num == 0) return "";
        String out = array[0].toString();
        for(int ii=1; ii < num; ii++) {
            out += (separator + array[ii].toString());
        }
        return out;
    }


    /**
     * This must exist in the jdk ?
     *
     * @param string
     * @return
     */
    public static int countChar(String string, char c) {
        int cnt = 0;
        for (int i = 0; i < string.length(); i++) {
            if (c == string.charAt(i)) {
                cnt++;
            }
        }
        return cnt;

    }

    /**
     * Strips quotation marks from string.  Quotes can be single or double, but must match each other.
     *
     * @param fileString
     * @return
     */
    public static String stripQuotes(String fileString) {
        if ((fileString.startsWith("\"") && fileString.endsWith("\"")) ||
                (fileString.startsWith("'") && fileString.endsWith("'"))) {
            fileString = fileString.substring(1, fileString.length() - 1);
        }
        return fileString;
    }

}
