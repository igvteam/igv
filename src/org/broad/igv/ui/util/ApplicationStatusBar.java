/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */


package org.broad.igv.ui.util;

import com.jidesoft.swing.JideBoxLayout;
import com.jidesoft.swing.JideButton;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.ui.FontManager;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;
import java.util.TimerTask;


/**
 * @author eflakes
 */
public class ApplicationStatusBar extends JPanel { //StatusBar {

    static Logger log = Logger.getLogger(ApplicationStatusBar.class);
    private JLabel messageBox;
    private JLabel messageBox2;
    private JLabel messageBox3;
    private JLabel memoryStatus;

    private JButton cancelButton;

    java.util.Timer timer;

    public ApplicationStatusBar() {
        initialize();
    }

    private void initialize() {

        setBackground(new Color(240, 240, 240));
        Color messageBG = new Color(230, 230, 230);
        Font messageFont = FontManager.getFont(11);

        setMinimumSize(new Dimension(200, 20));
        setPreferredSize(new Dimension(800, 20));

        JideBoxLayout layout = new JideBoxLayout(this, JideBoxLayout.X_AXIS);
        layout.setGap(3);
        setLayout(layout);

        messageBox = createMessageField(messageBG, messageFont);
        messageBox.setMinimumSize(new Dimension(135, 10));
        messageBox.setPreferredSize(new Dimension(135, 20));
        add(messageBox, JideBoxLayout.FIX);

        if (!Globals.isProduction()) {
            cancelButton = new JideButton(IconFactory.getInstance().getIcon(IconFactory.IconID.CLOSE));
            cancelButton.setMinimumSize(new Dimension(20, 10));
            cancelButton.setPreferredSize(new Dimension(20, 20));
            add(cancelButton, JideBoxLayout.FIX);
        }

        messageBox2 = createMessageField(messageBG, messageFont);
        messageBox2.setMinimumSize(new Dimension(150, 10));
        messageBox2.setPreferredSize(new Dimension(150, 20));
        add(messageBox2, JideBoxLayout.FIX);

        messageBox3 = createMessageField(messageBG, messageFont);
        messageBox3.setMinimumSize(new Dimension(165, 10));
        messageBox3.setPreferredSize(new Dimension(165, 20));
        add(messageBox3, JideBoxLayout.VARY);

        memoryStatus = createMessageField(messageBG, messageFont);
        memoryStatus.setPreferredSize(new Dimension(100, 20));
        memoryStatus.setMinimumSize(new Dimension(100, 10));
        memoryStatus.setBackground(messageBG);
        add(memoryStatus, JideBoxLayout.FIX);

        MemoryUpdateTask updateTask = new MemoryUpdateTask(memoryStatus);
        timer = new java.util.Timer();
        timer.schedule(updateTask, 0, 1000);


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


    private JLabel createMessageField(Color bg, Font font) {
        JLabel messageField = new JLabel();
        messageField.setBackground(bg);
        messageField.setFont(font);
        messageField.setBorder(BorderFactory.createLineBorder(Color.black));
        return messageField;

    }

    /**
     * Set the cancel button
     */
    public void activateCancelButton(ActionListener listener) {

        if (!Globals.isProduction()) {
            this.cancelButton.addActionListener(listener);
            this.cancelButton.addActionListener(new CancelButtonActionListener());
            this.cancelButton.setEnabled(true);
        }

    }

    public void deactivateCancelButton() {
        if (!Globals.isProduction()) {
            for (ActionListener l : this.cancelButton.getActionListeners()) {
                this.cancelButton.removeActionListener(l);
            }
            this.cancelButton.setEnabled(false);
        }
    }

    public JButton getCancelButton() {
        return cancelButton;
    }

    class CancelButtonActionListener implements ActionListener {
        @Override
        public void actionPerformed(ActionEvent e) {
            deactivateCancelButton();
        }
    }


    class MemoryUpdateTask extends TimerTask {

        JLabel textField;
        NumberFormat format;

        public MemoryUpdateTask(JLabel textField) {
            this.textField = textField;
            format = NumberFormat.getIntegerInstance();
        }

        @Override
        public void run() {
            Runtime runtime = Runtime.getRuntime();
            int freeMemory = (int) (runtime.freeMemory() / 1000000);
            int totalMemory = (int) (runtime.totalMemory() / 1000000);
            int usedMemory = (totalMemory - freeMemory);
            String um = format.format(usedMemory);
            String tm = format.format(totalMemory);
            textField.setText(um + "M of " + tm + "M");
        }

    }
}