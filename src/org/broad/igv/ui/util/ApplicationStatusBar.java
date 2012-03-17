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

import com.jidesoft.swing.JideBoxLayout;
import org.apache.log4j.Logger;
import org.broad.igv.ui.FontManager;

import javax.swing.*;
import java.text.NumberFormat;
import java.awt.*;
import java.util.TimerTask;


/**
 * @author eflakes
 */
public class ApplicationStatusBar extends JPanel { //StatusBar {

    static Logger log = Logger.getLogger(ApplicationStatusBar.class);
    private JTextField messageBox;
    private JTextField messageBox2;
    private JTextField messageBox3;
    private MemoryStatus memoryStatus;
    private Font font;
    //private LabelStatusBarItem messageBox3;
    //private MemoryStatusBarItem memoryBox;

    public ApplicationStatusBar() {
        initialize();
    }

    private void initialize() {

        Color bg = new Color(250, 250, 250);
        font = FontManager.getFont(10);

        setMinimumSize(new Dimension(200, 25));

        // setPreferredSize(new Dimension(800, 32));

        JideBoxLayout layout = new JideBoxLayout(this, JideBoxLayout.X_AXIS);
        layout.setGap(1);
        setLayout(layout);

        messageBox = new JTextField();
        messageBox.setMinimumSize(new Dimension(165, 10));
        messageBox.setPreferredSize(new Dimension(165, 20));
        messageBox.setEditable(false);
        messageBox.setFont(font);
        messageBox.setBackground(bg);
        add(messageBox, JideBoxLayout.FIX);

        messageBox2 = new JTextField();
        messageBox2.setMinimumSize(new Dimension(150, 10));
        messageBox2.setPreferredSize(new Dimension(150, 20));
        messageBox2.setEditable(false);
        messageBox2.setBackground(bg);
        messageBox2.setFont(font);
        add(messageBox2, JideBoxLayout.FIX);

        messageBox3 = new JTextField();
        messageBox3.setMinimumSize(new Dimension(165, 10));
        messageBox3.setPreferredSize(new Dimension(165, 20));
        messageBox3.setEditable(false);
        messageBox3.setBackground(bg);
        messageBox3.setFont(font);
        add(messageBox3, JideBoxLayout.VARY);

        memoryStatus = new MemoryStatus();
        memoryStatus.setPreferredSize(new Dimension(100, 20));
        memoryStatus.setMinimumSize(new Dimension(100, 10));
        memoryStatus.setBackground(bg);
        add(memoryStatus, JideBoxLayout.FIX);

    }

    public void setMessage(final String message) {
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                messageBox.setText(message);
                messageBox.paintImmediately(messageBox.getBounds());
            }
        });
    }

    public void setMessage2(final String message) {
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                messageBox2.setText(message);
                messageBox2.paintImmediately(messageBox2.getBounds());
            }
        });
    }


    class MemoryStatus extends JTextField {

        NumberFormat format;
        java.util.Timer timer;

        public MemoryStatus() {

            setFont(font);
            setEditable(false);
            format = NumberFormat.getIntegerInstance();

            TimerTask updateTask = new TimerTask() {
                @Override
                public void run() {
                    Runtime runtime = Runtime.getRuntime();
                    int freeMemory = (int) (runtime.freeMemory() / 1000000);
                    int totalMemory = (int) (runtime.totalMemory() / 1000000);
                    int usedMemory = (totalMemory - freeMemory);
                    String um = format.format(usedMemory);
                    String tm = format.format(totalMemory);
                    setText(um + "M of " + tm + "M");
                }
            };

            timer = new java.util.Timer();
            timer.schedule(updateTask, 0, 1000);


        }

    }
}
