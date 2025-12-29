/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui;

import org.igv.logging.*;

import javax.swing.*;
import java.awt.*;
import java.lang.Thread.UncaughtExceptionHandler;
import java.util.ConcurrentModificationException;

public class DefaultExceptionHandler implements UncaughtExceptionHandler {

    Logger log = LogManager.getLogger(DefaultExceptionHandler.class);

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
