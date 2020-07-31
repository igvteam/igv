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


package org.broad.igv.ui.util;

import com.jidesoft.swing.JideBoxLayout;
import com.jidesoft.swing.JideButton;
import org.apache.log4j.Logger;
import org.broad.igv.google.OAuthProvider;
import org.broad.igv.ui.FontManager;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.google.OAuthUtils;
//import org.broad.igv.event.StopEvent;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;
import java.util.TimerTask;


/**
 * @author eflakes
 */
public class ApplicationStatusBar extends JPanel implements IGVEventObserver { //StatusBar {

    static Logger log = Logger.getLogger(ApplicationStatusBar.class);
    public JButton stopButton;
    private JLabel messageBox;
    private JLabel messageBox2;
    private JLabel messageBox3;
    private JLabel memoryStatus;

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

        stopButton = new JButton();
        stopButton.setBorder(BorderFactory.createLineBorder(Color.black));
        stopButton.setIcon(new javax.swing.ImageIcon(
                getClass().getResource("/toolbarButtonGraphics/general/Stop16.gif")));
        stopButton.setMaximumSize(new java.awt.Dimension(16, 16));
        stopButton.setMinimumSize(new java.awt.Dimension(16, 16));
        stopButton.setPreferredSize(new java.awt.Dimension(16, 16));
        stopButton.setSize(new java.awt.Dimension(16, 16));
        stopButton.setToolTipText("Cancel load");
        stopButton.setEnabled(false);
        add(stopButton, JideBoxLayout.FIX);

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

        // Once we get positive result from the auth flow, perhaps fetch user full name and/or profile avatar?
        IGVEventBus.getInstance().subscribe(OAuthProvider.AuthStateEvent.class, this);
    }

    public void setMessage(final String message) {
        UIUtilities.invokeOnEventThread(() -> {
            messageBox.setText(message);
            messageBox.validate();
            messageBox.paintImmediately(messageBox.getBounds());
        });
    }

    public void setMessage2(final String message) {
        UIUtilities.invokeOnEventThread(() -> {
            messageBox2.setText(message);
            messageBox.validate();
            messageBox2.paintImmediately(messageBox2.getBounds());
        });
    }

    public void setMessage3(final String message) {
        UIUtilities.invokeOnEventThread(() -> {
            messageBox3.setText(message);
            messageBox.validate();
            messageBox3.paintImmediately(messageBox2.getBounds());
        });
    }

    private JLabel createMessageField(Color bg, Font font) {
        JLabel messageField = new JLabel("");
        messageField.setBackground(bg);
        messageField.setFont(font);
        messageField.setBorder(BorderFactory.createLineBorder(Color.black));
        return messageField;

    }

    public void enableStopButton(boolean enable) {
        stopButton.setEnabled(enable);
    }

    @Override
    public void receiveEvent(Object event) {
        if(event instanceof OAuthProvider.AuthStateEvent) {
            boolean isAuthed = ((OAuthProvider.AuthStateEvent) event).isAuthenticated();
            String authProvider = ((OAuthProvider.AuthStateEvent) event).getAuthProvider();
            String userName = ((OAuthProvider.AuthStateEvent) event).getUserName();
            String email = ((OAuthProvider.AuthStateEvent) event).getEmail();

            if (isAuthed) {
                setMessage3("Logged in as: " + email + " via " + authProvider);
            } else {
                setMessage3("");
            }
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
            UIUtilities.invokeOnEventThread(() -> textField.setText(um + "M of " + tm + "M"));
        }

    }
}
