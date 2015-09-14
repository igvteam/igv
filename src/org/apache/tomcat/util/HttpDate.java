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

import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.text.*;

/**
 * This class can be used to efficiently parse and write an RFC 1123
 * formatted date in an HTTP message header.  Also supports reading the
 * RFC 1036 format and ANSI C's asctime() format, as suggested by HTTP/1.0
 * and mandated by HTTP/1.1.
 *
 * @author dac@eng.sun.com
 * @author Jason Hunter [jch@eng.sun.com]
 * @author James Todd [gonzo@eng.sun.com]
 */

public class HttpDate extends Ascii {

//    private StringManager sm =
//        StringManager.getManager(Constants.Package);

    // ONLY FOR COMPAT -- KILL ASAP -- just make sure that dependant
    // classes know what's up. ref. MimeHeaderField
    private static final String DATESTR = "Sun, 06 Nov 1994 08:49:37 GMT";
    public static final int DATELEN = DATESTR.length();
    // END COMPAT -- DON'T FORGET TO KILL

    // we force our locale here as all http dates are in english
    private final static Locale loc = Locale.US;

    // all http dates are expressed as time at GMT
    private final static TimeZone zone = TimeZone.getTimeZone("GMT");

    // format for RFC 1123 date string -- "Sun, 06 Nov 1994 08:49:37 GMT"
    private final static String rfc1123Pattern =
            "EEE, dd MMM yyyy HH:mm:ss z";

    // format for RFC 1036 date string -- "Sunday, 06-Nov-94 08:49:37 GMT"
    private final static String rfc1036Pattern =
            "EEEEEEEEE, dd-MMM-yy HH:mm:ss z";

    // format for C asctime() date string -- "Sun Nov  6 08:49:37 1994"
    private final static String asctimePattern =
            "EEE MMM d HH:mm:ss yyyy";

    private final static SimpleDateFormat rfc1123Format =
            new SimpleDateFormat(rfc1123Pattern, loc);

    private final static SimpleDateFormat rfc1036Format =
            new SimpleDateFormat(rfc1036Pattern, loc);

    private final static SimpleDateFormat asctimeFormat =
            new SimpleDateFormat(asctimePattern, loc);

    static {
        rfc1123Format.setTimeZone(zone);
        rfc1036Format.setTimeZone(zone);
        asctimeFormat.setTimeZone(zone);
    }

    // protected so that oldcookieexpiry in cookieutils can use
    // yes, this is sloppy as crap and could stand to be done better.
    protected Calendar calendar = new GregorianCalendar(zone, loc);

    public HttpDate() {
        calendar.setTime(new Date(System.currentTimeMillis()));
    }

    public HttpDate(long ms) {
        calendar.setTime(new Date(ms));
    }

    public void setTime() {
        calendar.setTime(new Date(System.currentTimeMillis()));
    }

    public void setTime(long ms) {
        calendar.setTime(new Date(ms));
    }

    public void parse(String dateString) {
        try {
            Date date = rfc1123Format.parse(dateString);
            calendar.setTime(date);
            return;
        } catch (ParseException e) {
        }

        try {
            Date date = rfc1036Format.parse(dateString);
            calendar.setTime(date);
            return;
        } catch (ParseException e) {
        }

        try {
            Date date = asctimeFormat.parse(dateString);
            calendar.setTime(date);
            return;
        } catch (ParseException pe) {
            //String msg = sm.getString("httpDate.pe", dateString);
            String msg = "Could not parse data: " + dateString;
            throw new IllegalArgumentException(msg);
        }
    }

    public void parse(byte[] b, int off, int len) {
        // ok -- so this is pretty stoopid, but the old version of this
        // source took this arg set, so we will too for now (backwards compat)
        String dateString = new String(b, off, len);
        parse(dateString);
    }

    public void write(OutputStream out) throws IOException {
        String dateString = rfc1123Format.format(calendar.getTime());
        byte[] b = dateString.getBytes();
        out.write(b);
    }

    public String toString() {
        return rfc1123Format.format(calendar.getTime());
    }

    public long getTime() {
        return calendar.getTime().getTime();
    }

    public static long getCurrentTime() {
        return System.currentTimeMillis();
    }
//
//    // KILL, THIS IS ONLY HERE FOR TEMP COMPAT as MimeHeaderField uses it.
//    public int getBytes(byte[] buf, int off, int len) {
//	if (len < DATELEN) {
//            String msg = sm.getString("httpDate.iae", new Integer(len));
//
//	    throw new IllegalArgumentException(msg);
//	}
//
//	String dateString = rfc1123Format.format(calendar.getTime());
//	byte[] b = dateString.getBytes();
//	System.arraycopy(b, 0, buf, off, DATELEN);
//	return DATELEN;
//    }
}
