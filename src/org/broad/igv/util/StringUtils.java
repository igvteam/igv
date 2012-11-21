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

package org.broad.igv.util;

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
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

        ArrayList<String> strings = new ArrayList();
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
     * support UTF-8 (it's part of the Java spec), URLDecoder.decode
     * is used.
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

    public static String encodeURL(String s) {
        if (s == null) {
            return null;
        }
        try {
            return URLEncoder.encode(s, "UTF-8");
        } catch (UnsupportedEncodingException e) {
            throw new RuntimeException(e);
        }
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
}
