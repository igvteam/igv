/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.util;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;

/**
 * @author jrobinso
 */
public class MessageUtils {

    private static Logger log = Logger.getLogger(MessageUtils.class);

    public static void showMessage(String message) {

        if (Globals.isHeadless() || Globals.isSuppressMessages() || !IGV.hasInstance()) {
            log.info(message);
        } else {
            JOptionPane.showMessageDialog(IGV.getMainFrame(), message);
        }
    }


    public static boolean confirm(String message) {

        return confirm(IGV.getMainFrame(), message);

    }

    public static boolean confirm(Component component, String message) {

        int opt = JOptionPane.showConfirmDialog(component, message, "Confirm", JOptionPane.YES_NO_OPTION);
        return opt == JOptionPane.YES_OPTION;

    }

    /**
     * Method description
     *
     * @param component
     * @param message
     * @param log
     * @param e
     */
    public static void showAndLogErrorMessage(final Component component, final String message,
                                              final Logger log, final Exception e) {
        if (log != null) {

            if (e != null) {
                log.error(message, e);
            } else {
                log.error(message);
            }
        }
        JOptionPane.showMessageDialog(component, message);

    }

    /**
     * Method description
     *
     * @param component
     * @param message
     * @param log
     */
    public static void showAndLogErrorMessage(Component component, String message, Logger log) {

        showAndLogErrorMessage(component, message, log, null);
    }

    public static String showInputDialog(String message, String defaultValue) {
        String val = JOptionPane.showInputDialog(IGV.getMainFrame(), message, defaultValue);
        return val;
    }

    public static String showInputDialog(String message) {
        String val = JOptionPane.showInputDialog(IGV.getMainFrame(), message);
        return val;
    }
}
