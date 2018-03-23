/*
*
* The Apache Software License, Version 1.1
*
* Copyright (c) 1999 The Apache Software Foundation.  All rights
* reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in
*    the documentation and/or other materials provided with the
*    distribution.
*
* 3. The end-user documentation included with the redistribution, if
*    any, must include the following acknowlegement:
*       "This product includes software developed by the
*        Apache Software Foundation (http://www.apache.org/)."
*    Alternately, this acknowlegement may appear in the software itself,
*    if and wherever such third-party acknowlegements normally appear.
*
* 4. The names "The Jakarta Project", "Tomcat", and "Apache Software
*    Foundation" must not be used to endorse or promote products derived
*    from this software without prior written permission. For written
*    permission, please contact apache@apache.org.
*
* 5. Products derived from this software may not be called "Apache"
*    nor may "Apache" appear in their names without prior written
*    permission of the Apache Group.
*
* THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
* OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED.  IN NO EVENT SHALL THE APACHE SOFTWARE FOUNDATION OR
* ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
* OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
* SUCH DAMAGE.
* ====================================================================
*
* This software consists of voluntary contributions made by many
* individuals on behalf of the Apache Software Foundation.  For more
* information on the Apache Software Foundation, please see
* <http://www.apache.org/>.
*
* [Additional notices, if required by prior licensing conditions]
*
*/


package org.apache.tomcat.util;

/**
 * This class implements some basic ASCII character handling functions.
 *
 * @author dac@eng.sun.com
 * @author James Todd [gonzo@eng.sun.com]
 */

public class Ascii {
    /*
     * Character translation tables.
     */

    private static final byte[] toUpper = new byte[256];
    private static final byte[] toLower = new byte[256];

    /*
     * Character type tables.
     */

    private static final boolean[] isAlpha = new boolean[256];
    private static final boolean[] isUpper = new boolean[256];
    private static final boolean[] isLower = new boolean[256];
    private static final boolean[] isWhite = new boolean[256];
    private static final boolean[] isDigit = new boolean[256];

    /*
     * Initialize character translation and type tables.
     */

    static {
        for (int i = 0; i < 256; i++) {
            toUpper[i] = (byte) i;
            toLower[i] = (byte) i;
        }

        for (int lc = 'a'; lc <= 'z'; lc++) {
            int uc = lc + 'A' - 'a';

            toUpper[lc] = (byte) uc;
            toLower[uc] = (byte) lc;
            isAlpha[lc] = true;
            isAlpha[uc] = true;
            isLower[lc] = true;
            isUpper[uc] = true;
        }

        isWhite[' '] = true;
        isWhite['\t'] = true;
        isWhite['\r'] = true;
        isWhite['\n'] = true;
        isWhite['\f'] = true;
        isWhite['\b'] = true;

        for (int d = '0'; d <= '9'; d++) {
            isDigit[d] = true;
        }
    }

    /**
     * Returns the upper case equivalent of the specified ASCII character.
     */

    public static int toUpper(int c) {
        return toUpper[c & 0xff] & 0xff;
    }

    /**
     * Returns the lower case equivalent of the specified ASCII character.
     */

    public static int toLower(int c) {
        return toLower[c & 0xff] & 0xff;
    }

    /**
     * Returns true if the specified ASCII character is upper or lower case.
     */

    public static boolean isAlpha(int c) {
        return isAlpha[c & 0xff];
    }

    /**
     * Returns true if the specified ASCII character is upper case.
     */

    public static boolean isUpper(int c) {
        return isUpper[c & 0xff];
    }

    /**
     * Returns true if the specified ASCII character is lower case.
     */

    public static boolean isLower(int c) {
        return isLower[c & 0xff];
    }

    /**
     * Returns true if the specified ASCII character is white space.
     */

    public static boolean isWhite(int c) {
        return isWhite[c & 0xff];
    }

    /**
     * Returns true if the specified ASCII character is a digit.
     */

    public static boolean isDigit(int c) {
        return isDigit[c & 0xff];
    }

    /**
     * Parses an unsigned integer from the specified subarray of bytes.
     *
     * @param b   the bytes to parse
     * @param off the start offset of the bytes
     * @param len the length of the bytes
     * @throws NumberFormatException if the integer format was invalid
     */

    public static int parseInt(byte[] b, int off, int len)
            throws NumberFormatException {
        int c;

        if (b == null || len <= 0 || !isDigit(c = b[off++])) {
//            StringManager sm =
//                StringManager.getManager(Constants.Package);
//            String msg = sm.getString("ascii.parseInit.nfe", b);
            String msg = "Could not parse byte array as integer: " + new String(b);
            throw new NumberFormatException(msg);
        }

        int n = c - '0';

        while (--len > 0) {
            if (!isDigit(c = b[off++])) {
//                StringManager sm =
//                    StringManager.getManager(Constants.Package);
                //    String msg = sm.getString("ascii.parseInit.nfe", b);
                String msg = "Could not parse byte array as integer: " + new String(b);

                throw new NumberFormatException(msg);
            }

            n = n * 10 + c - '0';
        }

        return n;
    }
}
