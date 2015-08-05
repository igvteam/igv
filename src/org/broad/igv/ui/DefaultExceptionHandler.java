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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.lang.Thread.UncaughtExceptionHandler;
import java.util.ConcurrentModificationException;

public class DefaultExceptionHandler implements UncaughtExceptionHandler {

    Logger log = Logger.getLogger(DefaultExceptionHandler.class);

    public void uncaughtException(Thread t, Throwable e) {
        if (e instanceof ConcurrentModificationException) {
            // Ignore these,  they are logged elsewhere
        } else {
            //JOptionPane.showMessageDialog(findActiveFrame(),
            //        "An unexpected error occured: " + e.toString(), "Exception Occurred", JOptionPane.OK_OPTION);
            log.error("Unhandled exception", e);

        }
    }

    private Frame findActiveFrame() {
        Frame[] frames = JFrame.getFrames();
        for (int i = 0; i < frames.length; i++) {
            if (frames[i].isVisible()) {
                return frames[i];
            }
        }
        return null;
    }
}
