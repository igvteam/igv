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
package org.broad.igv.ui.util;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.lang.reflect.InvocationTargetException;

/**
 * Provides thread-safe, Swing-safe, utilities for interacting with JOptionPane.  Accounts for
 * (1) Swing is not thread safe => synchronize access
 * (2) JOptionPane methods must be invoked on event dispatch thread
 *
 * @author jrobinso
 */
public class MessageUtils {

    private static Logger log = Logger.getLogger(MessageUtils.class);

    // Somewhat silly class, needed to pass values between threads
    static class ValueHolder {
        Object value;
    }

    /**
     * Log the exception and show {@code message} to the user
     *
     * @param e
     * @param message
     */
    public static void showErrorMessage(String message, Exception e) {
        log.error(message, e);
        showMessage(Level.ERROR, message);
    }

    public static void showMessage(String message) {
        showMessage(Level.INFO, message);
    }

    public static synchronized void showMessage(Level level, String message) {

        log.log(level, message);
        boolean showDialog = !(Globals.isHeadless() || Globals.isSuppressMessages() || Globals.isTesting());
        if (showDialog) {
            // Always use HTML for message displays, but first remove any embedded <html> tags.
            message = "<html>" + message.replaceAll("<html>", "");
            Frame parent = IGV.hasInstance() ? IGV.getMainFrame() : null;
            Color background = parent != null ? parent.getBackground() : Color.lightGray;

            //JEditorPane So users can select text
            JEditorPane content = new JEditorPane();
            content.setContentType("text/html");
            content.setText(message);
            content.setBackground(background);
            content.setEditable(false);
            Component dispMessage = content;

            //Really long messages should be scrollable
            if(message.length() > 200){
                Dimension size = new Dimension(1000, content.getHeight());
                content.setPreferredSize(size);
                JScrollPane pane = new JScrollPane(content);
                dispMessage = pane;
            }

            JOptionPane.showMessageDialog(parent, dispMessage);
        }
    }

    public static void setStatusBarMessage(final String message) {
        log.debug("Status bar: " + message);
        if (IGV.hasInstance()) {
            IGV.getInstance().setStatusBarMessage(message);
        }
    }

    public static synchronized boolean confirm(final String message) {
        if(Globals.isHeadless()){
            log.error("Attempted to confirm while running headless with the following message:\n" + message);
            return false;
        }

        final Frame parent = IGV.hasInstance() ? IGV.getMainFrame() : null;
        return confirm(parent, message);
    }

    /**
     * Show a yes/no confirmation dialog.
     *
     * @param component
     * @param message
     * @return
     */
    public static synchronized boolean confirm(final Component component, final String message) {


        if (SwingUtilities.isEventDispatchThread()) {
            int opt = JOptionPane.showConfirmDialog(component, message, "Confirm", JOptionPane.YES_NO_OPTION);
            return opt == JOptionPane.YES_OPTION;
        } else {
            final ValueHolder returnValue = new ValueHolder();
            Runnable runnable = new Runnable() {
                public void run() {
                    int opt = JOptionPane.showConfirmDialog(component, message, "Confirm", JOptionPane.YES_NO_OPTION);
                    returnValue.value = (opt == JOptionPane.YES_OPTION);
                }
            };
            try {
                SwingUtilities.invokeAndWait(runnable);
            } catch (InterruptedException e) {
                log.error("Error in showMessage", e);
                throw new RuntimeException(e);
            } catch (InvocationTargetException e) {
                log.error("Error in showMessage", e);
                throw new RuntimeException(e.getCause());
            }

            return (Boolean) (returnValue.value);

        }
    }

    public static String showInputDialog(String message, final String defaultValue) {

        final Frame parent = IGV.hasInstance() ? IGV.getMainFrame() : null;
        //Pad message with spaces so it's as wide as the defaultValue
        if(message.length() < defaultValue.length()){
            message = String.format("%-" + defaultValue.length() + "s", message);
        }
        final String actMsg = message;

        if (SwingUtilities.isEventDispatchThread()) {
            String val = JOptionPane.showInputDialog(parent, actMsg, defaultValue);
            return val;
        } else {
            final ValueHolder returnValue = new ValueHolder();
            Runnable runnable = new Runnable() {
                public void run() {
                    String val = JOptionPane.showInputDialog(parent, actMsg, defaultValue);
                    returnValue.value = val;
                }
            };
            try {
                SwingUtilities.invokeAndWait(runnable);
            } catch (InterruptedException e) {
                log.error("Error in showMessage", e);
                throw new RuntimeException(e);
            } catch (InvocationTargetException e) {
                log.error("Error in showMessage", e);
                throw new RuntimeException(e.getCause());
            }

            return (String) (returnValue.value);
        }
    }

    public static String showInputDialog(final String message) {

        final Frame parent = IGV.hasInstance() ? IGV.getMainFrame() : null;
        if (SwingUtilities.isEventDispatchThread()) {
            String val = JOptionPane.showInputDialog(parent, message);
            return val;
        } else {
            final ValueHolder returnValue = new ValueHolder();
            Runnable runnable = new Runnable() {
                public void run() {
                    String val = JOptionPane.showInputDialog(parent, message);
                    returnValue.value = val;
                }
            };
            try {
                SwingUtilities.invokeAndWait(runnable);
            } catch (InterruptedException e) {
                log.error("Error in showMessage", e);
                throw new RuntimeException(e);
            } catch (InvocationTargetException e) {
                log.error("Error in showMessage", e);
                throw new RuntimeException(e.getCause());
            }

            return (String) (returnValue.value);
        }
    }


    /**
     * Test program - call all methods from both main and swing threads
     *
     * @param args
     * @throws Exception
     */

    public static void main(String[] args) throws Exception {

        Runnable runnable = new Runnable() {
            public void run() {
                showMessage("showMessage");

                confirm("confirm");

                confirm(null, "confirm with parent");

                showInputDialog("showInputDialog", "default");

                showInputDialog("showInputDialog");
            }
        };

        // Test on main thread
        runnable.run();


        // Test on swing thread
        SwingUtilities.invokeLater(runnable);

    }

}
