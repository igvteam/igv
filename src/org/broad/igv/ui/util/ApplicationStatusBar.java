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


package org.broad.igv.ui.util;

//~--- non-JDK imports --------------------------------------------------------


//import com.jidesoft.status.LabelStatusBarItem;
//import com.jidesoft.status.MemoryStatusBarItem;
//import com.jidesoft.status.StatusBar;
//import com.jidesoft.swing.JideBoxLayout;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.text.DecimalFormat;


/**
 * @author eflakes
 */
public class ApplicationStatusBar extends JPanel { //StatusBar {

    static Logger log = Logger.getLogger(ApplicationStatusBar.class);
    //private LabelStatusBarItem messageBox;
    //private LabelStatusBarItem messageBox2;
    //private LabelStatusBarItem messageBox3;
    //private MemoryStatusBarItem memoryBox;

    public ApplicationStatusBar() {
        initialize();
    }

    private void initialize() {

//        messageBox = new LabelStatusBarItem("Line");
//        messageBox.setAlignment(JLabel.LEFT);
//        messageBox.setMinimumSize(new Dimension(165, 10));
//        add(messageBox, JideBoxLayout.FLEXIBLE);
//
//        messageBox2 = new LabelStatusBarItem("Line");
//        messageBox2.setAlignment(JLabel.LEFT);
//        messageBox2.setMinimumSize(new Dimension(165, 10));
//        add(messageBox2, JideBoxLayout.FLEXIBLE);
//
//        messageBox3 = new LabelStatusBarItem("Line");
//        messageBox3.setHorizontalAlignment(SwingConstants.RIGHT);
//        add(messageBox3, JideBoxLayout.VARY);
//
//        memoryBox = new MemoryStatusBarItem();
//
//        add(memoryBox, JideBoxLayout.FLEXIBLE);
    }

    public void setMessage(final String message) {
 //       UIUtilities.invokeOnEventThread(new Runnable() {
 //
//            public void run() {
//                messageBox.setText(message);
//                messageBox.paintImmediately(messageBox.getBounds());
//            }
 //       });
    }

    public void setMessage2(final String message) {
//        UIUtilities.invokeOnEventThread(new Runnable() {
//
//            public void run() {
//                messageBox2.setText(message);
//                messageBox2.paintImmediately(messageBox2.getBounds());
//            }
//        });
    }
}
