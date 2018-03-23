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

import org.apache.log4j.Logger;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;

/**
 * @author Joshua Gould
 */
public class BrowserLauncher {

    static Logger log = Logger.getLogger(BrowserLauncher.class);

    /**
     * The Java virtual machine that we are running on. Actually, in most cases we only care about the operating system,
     * but some operating systems require us to switch on the VM.
     */
    private static int jvm;

    /**
     * JVM constant for any Windows JVM
     */
    private static final int WINDOWS = 1;

    /**
     * JVM constant for any Mac JVM
     */
    private static final int MAC = 2;

    /**
     * JVM constant for any other platform
     */
    private static final int OTHER = 3;

    /**
     * specified the path to the browser
     */
    private static final int SPECIFIED = 4;

    /**
     * The browser for the system
     */
    private static String specifiedBrowser;

    /**
     * An initialization block that determines the operating system.
     */
    static {
        String osName = System.getProperty("os.name");
        if (osName.startsWith("Mac OS")) {
            jvm = MAC;
        } else if (osName.startsWith("Windows")) {
            jvm = WINDOWS;
        } else {
            jvm = OTHER;
        }
    }

    /**
     * This class should be never be instantiated; this just ensures so.
     */
    private BrowserLauncher() {
    }

    /**
     * Attempts to open the default web browser to the given URL.
     *
     * @param url The URL to open
     * @throws IOException If the web browser could not be located or does not run
     */
    public static void openURL(String url) throws IOException {
        if (!(Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE))) {
            openURL_old(url);
        } else {
            try {
                Desktop.getDesktop().browse(new URI(url));
            } catch (URISyntaxException e) {
                log.error("Error opening url " + url, e);
            }

        }
    }

    /**
     * Attempts to open the default web browser to the given URL.
     *
     * @param url The URL to open
     * @throws IOException If the web browser could not be located or does not run
     */
    public static void openURL(URL url) throws IOException {
        if (!(Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE))) {
            openURL_old(url.toExternalForm());
        } else {
            try {
                Desktop.getDesktop().browse(url.toURI());
            } catch (URISyntaxException e) {
                log.error("Error opening url " + url, e);
            }

        }
    }

    private static void openURL_old(String url) throws IOException {
        switch (jvm) {
            case MAC:
                Runtime.getRuntime().exec(new String[]{"/usr/bin/open", url});
                break;

            case WINDOWS:
                Process process = Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + url);
                // This avoids a memory leak on some versions of Java on Windows.
                // That's hinted at in
                // <http://developer.java.sun.com/developer/qow/archive/68/>.
                try {
                    process.waitFor();
                    process.exitValue();
                } catch (InterruptedException ie) {
                    throw new IOException("InterruptedException while launching browser: " + ie.getMessage());
                }
                break;
            case OTHER:
                if (new File("/usr/bin/gnome-open").exists()) {
                    Runtime.getRuntime().exec(new String[]{"/usr/bin/gnome-open", url});
                } else if (new File("/usr/bin/kde-open").exists()) {
                    Runtime.getRuntime().exec(new String[]{"/usr/bin/kde-open", url});
                }

                break;
            case SPECIFIED:
                process = Runtime.getRuntime().exec(new String[]{specifiedBrowser, url});
                try {
                    process.waitFor();
                    process.exitValue();
                } catch (InterruptedException ie) {
                    throw new IOException("InterruptedException while launching browser: " + ie.getMessage());
                }
                break;
        }
    }

    public static void setSpecifiedBrowser(String s) {
        specifiedBrowser = s;
        jvm = SPECIFIED;
    }

}