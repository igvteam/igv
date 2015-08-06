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

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.lang.reflect.InvocationTargetException;


/**
 * @author eflakes
 */
public class UIUtilities {

    final private static Logger log = Logger.getLogger(UIUtilities.class);

    final private static StringBuffer scratchBuffer = new StringBuffer();

    /**
     * Display a dialog which can be used to select a color
     *
     * @param dialogTitle
     * @param defaultColor The currently selected color
     * @return  The color the user selected, or null if none/cancelled
     */
    public static Color showColorChooserDialog(String dialogTitle, Color defaultColor) {

        Color color = null;
        JColorChooser chooser = new JColorChooser();
        chooser.setColor(defaultColor);
        while (true) {

            int response = JOptionPane.showConfirmDialog(IGV.getMainFrame(), chooser,
                    dialogTitle, JOptionPane.OK_CANCEL_OPTION);

            if ((response == JOptionPane.CANCEL_OPTION) || (response == JOptionPane.CLOSED_OPTION)) {
                return null;
            }

            color = chooser.getColor();
            if (color == null) {
                continue;
            } else {
                break;
            }
        }
        return color;
    }

    /**
     * Method description
     *
     * @param parent
     * @param message
     * @return
     */
    public static boolean showConfirmationDialog(Component parent, String message) {

        int status = JOptionPane.showConfirmDialog(parent, message, null,
                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE, null);

        if ((status == JOptionPane.CANCEL_OPTION) || (status == JOptionPane.CLOSED_OPTION)) {
            return false;
        }
        return true;
    }

    /**
     * Method description
     *
     * @param color
     * @return
     */
    public static String getcommaSeparatedRGBString(Color color) {

        if (color != null) {

            scratchBuffer.delete(0, scratchBuffer.length());    // Clear
            int red = color.getRed();
            int green = color.getGreen();
            int blue = color.getBlue();
            scratchBuffer.append(red);
            scratchBuffer.append(",");
            scratchBuffer.append(green);
            scratchBuffer.append(",");
            scratchBuffer.append(blue);
        }
        return scratchBuffer.toString();

    }

    /**
     * Method description
     *
     * @param window
     */
    public static void centerWindow(Window window) {

        Dimension dimension = window.getSize();
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        int x = (screenSize.width - dimension.width) / 2;
        int y = (screenSize.height - dimension.height) / 2;
        window.setLocation(x, y);
        window.requestFocus();
    }

    /**
     * A wrapper around invokeOnEventThread.  If the runnable is already in the event dispatching
     * queue it is just run.  Otherwise it is placed in the queue via invokeOnEventThread.
     * <p/>
     * I'm not sure this is strictly necessary,  but is safe.
     *
     * @param runnable
     */
    public static void invokeOnEventThread(Runnable runnable) {
        if (SwingUtilities.isEventDispatchThread()) {
            runnable.run();
        } else {
            SwingUtilities.invokeLater(runnable);
        }
    }

    /**
     * A wrapper around invokeOnEventThread.  If the runnable is already in the event dispatching
     * queue it is just run.  Otherwise it is placed in the queue via invokeOnEventThread.
     * <p/>
     * I'm not sure this is strictly necessary,  but is safe.
     *
     * @param runnable
     */
    public static void invokeAndWaitOnEventThread(Runnable runnable) {
        if (SwingUtilities.isEventDispatchThread()) {
            runnable.run();
        } else {
            try {
                SwingUtilities.invokeAndWait(runnable);
            } catch (InterruptedException e) {
                log.error("Error invoking runnable", e);
            } catch (InvocationTargetException e) {
                log.error("Error invoking runnable", e);
            }
        }
    }
}
