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

import java.awt.Desktop;
import java.awt.GraphicsEnvironment;
import java.awt.Image;
import java.awt.Taskbar;

import javax.swing.JOptionPane;

/**
 * Java version-specific integration with OS X (macOS)
 * @author eby
 */
public class DesktopIntegration {
    public static final void verifyJavaPlatform() {
        String javaVersion = System.getProperty("java.version");
        if (javaVersion == null || javaVersion.startsWith("1.8")) {
            try {
                System.out.println("Detected an unsupported Java version.  Java 8 is not supported by this release.");

                if (!GraphicsEnvironment.isHeadless()) {
                    JOptionPane.showMessageDialog(null, "Detected an unsupported Java version.  Java 8 is not supported by this release.");
                }
            } finally {
                System.exit(1);
            }
        }
    }
    
    public static void setDockIcon(Image image) {
        Taskbar.getTaskbar().setIconImage(image);
    }

    public static void setAboutHandler(IGVMenuBar igvMenuBar) {
        Desktop.getDesktop().setAboutHandler(e -> igvMenuBar.showAboutDialog());
    }
    
    public static void setQuitHandler() {
        Desktop.getDesktop().setQuitHandler((e, response) -> {
            try {
                ShutdownThread.runS();
            } finally {
                response.performQuit();
            }
        });
    }
}
