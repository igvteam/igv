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

import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.variant.SelectVcfFieldDialog;
import java.util.List;
import javax.swing.*;
import java.awt.*;
import java.lang.reflect.InvocationTargetException;
import java.util.Optional;
import java.util.function.Supplier;

/**
 * Provides thread-safe, Swing-safe, utilities for interacting with JOptionPane.  Accounts for
 * (1) Swing is not thread safe => synchronize access
 * (2) JOptionPane methods must be invoked on event dispatch thread
 *
 * @author jrobinso
 */
public class MessageUtils {

    private static Logger log = LogManager.getLogger(MessageUtils.class);

    // Somewhat silly class, needed to pass values between threads
    private static class ValueHolder<T> {
        T value;
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
        showMessage(null, message);
    }

    public static void showMessage(Level level, String message) {

        if(level != null) log.log(level, message);
        boolean showDialog = !(Globals.isHeadless() || Globals.isSuppressMessages() || Globals.isTesting() || Globals.isBatch());
        if (showDialog) {
            UIUtilities.invokeAndWaitOnEventThread(() -> {
                // Always use HTML for message displays, but first remove any embedded <html> tags.
                String dlgMessage = "<html>" + message.replaceAll("<html>", "");
                Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
                Color background = parent != null ? parent.getBackground() : Color.lightGray;

                //JEditorPane So users can select text
                JEditorPane content = new JEditorPane();
                content.setContentType("text/html");
                content.setText(dlgMessage);
                content.setBackground(background);
                content.setEditable(false);
                Component dispMessage = content;

                //Really long messages should be scrollable
                if (dlgMessage.length() > 200) {
                    Dimension size = new Dimension(1000, content.getHeight() + 100);
                    content.setPreferredSize(size);
                    JScrollPane pane = new JScrollPane(content);
                    dispMessage = pane;
                }

                JOptionPane.showMessageDialog(parent, dispMessage);
            });
        }
    }

    public static void setStatusBarMessage(final String message) {
        log.debug("Status bar: " + message);
        if (IGV.hasInstance()) {
            IGV.getInstance().setStatusBarMessage(message);
        }
    }

    public static synchronized boolean confirm(final String message) {
        if (Globals.isHeadless()) {
            log.error("Attempted to confirm while running headless with the following message:\n" + message);
            return true;
        }

        if (Globals.isBatch()) {
            return true;
        }

        final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
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

        if (Globals.isHeadless() || Globals.isBatch()) {
            return true;
        }

        return threadAwareInvokeAndGetValue(()  -> {
            int opt = JOptionPane.showConfirmDialog(component, message, "Confirm", JOptionPane.YES_NO_OPTION);
            return opt == JOptionPane.YES_OPTION;
        });
    }

    ;
    public static Optional<SelectVcfFieldDialog.ColorResult> showInputDialog(String message, String defaultValue, List<? extends VCFCompoundHeaderLine> valueOptions) {

        final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
        //Pad message with spaces so it's as wide as the defaultValue
        final String actMsg = padToMatchWidth(message, defaultValue);
        return threadAwareInvokeAndGetValue(() ->  SelectVcfFieldDialog.showValueChooseDialog(parent, actMsg, defaultValue, valueOptions));
    }

    public static String showInputDialog(String message, String defaultValue) {

        final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
        final String actMsg = padToMatchWidth(message, defaultValue);
        return threadAwareInvokeAndGetValue(() -> JOptionPane.showInputDialog(parent, actMsg, defaultValue));
    }

    private static String padToMatchWidth(String message, final String defaultValue) {
        //Pad message with spaces so it's as wide as the defaultValue
        if (defaultValue != null && message.length() < defaultValue.length()) {
            message = String.format("%-" + defaultValue.length() + "s", message);
        }
        return message;
    }

    /**
     * If this is the event thread run the supplier here, otherwise run it on the event thread
     */
    private static <T> T threadAwareInvokeAndGetValue(Supplier<T> supplier){
        if (SwingUtilities.isEventDispatchThread()) {
            return supplier.get();
        } else {
            return invokeAndGetValue(supplier);
        }
    }

    /**
     * Wraps a Supplier into a Runnable and pass it to SwingUtilities.invokeAndWait and then return the produced value
     * Useful for popping up dialogue boxes that provide a result.
     * @param supplier a function which must be run on the Swing Event Thread
     * @return the result of supplier.get(), throws a RuntimeException if an exception occurs
     * @param <T>  type of the returned value
     */
    private synchronized static <T> T invokeAndGetValue(Supplier<T> supplier) {
        final ValueHolder<T> returnValue = new ValueHolder<>();
        Runnable runnable = () -> returnValue.value = supplier.get();
        try {
            SwingUtilities.invokeAndWait(runnable);
        } catch (InterruptedException e) {
            log.error("Error in showMessage", e);
            throw new RuntimeException(e);
        } catch (InvocationTargetException e) {
            log.error("Error in showMessage", e);
            throw new RuntimeException(e.getCause());
        }

        return returnValue.value;
    }

    public static String showInputDialog(final String message) {

        final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
        return threadAwareInvokeAndGetValue(() -> JOptionPane.showInputDialog(parent, message));
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
