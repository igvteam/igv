/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2018 Broad Institute
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
package org.broad.igv.ui;

import java.awt.Image;

import org.apache.log4j.Logger;

import apple.dts.samplecode.osxadapter.OSXAdapter;

/**
 * Java version-specific integration with OS X (macOS)
 * @author eby
 */
public class OSXIntegration {
    private static Logger log = Logger.getLogger(OSXIntegration.class);
    
    public static void setDockIcon(Image image) {
        OSXAdapter.setDockIconImage(image);
    }

    public static void setAboutHandler(IGVMenuBar igvMenuBar) {
        try {
            OSXAdapter.setAboutHandler(igvMenuBar, igvMenuBar.getClass().getDeclaredMethod("showAboutDialog", (Class[]) null));
        } catch (Exception e) {
            log.error("Error setting apple-specific about handler", e);
        }
    }
    
    public static void setQuitHandler() {
        try {
            OSXAdapter.setQuitHandler(ShutdownThread.class, ShutdownThread.class.getDeclaredMethod("runS", (Class[]) null));
        } catch (Exception e) {
            log.error("Error setting apple-specific quit handler", e);
        }
    }
}
