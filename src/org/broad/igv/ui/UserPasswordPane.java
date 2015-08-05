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

package org.broad.igv.ui;

import org.broad.igv.util.StringUtils;

import javax.swing.*;
import java.awt.*;

/**
 * @author jrobinso
 * @date Mar 2, 2011
 */
public class UserPasswordPane {

    JPanel passPanel;
    JTextField userField;
    JPasswordField passwordField;
    String resourceString;

    public UserPasswordPane(String resourceString) {
        this.resourceString = resourceString;
        init();
    }

    private void init() {
        passPanel = new JPanel();
        passPanel.setLayout(new GridLayout(6, 1));

        JLabel message = new JLabel("Please enter your username and password for:");
        JLabel location = new JLabel(StringUtils.checkLength(resourceString, 80));
        JLabel username = new JLabel("User:");
        JLabel password = new JLabel("Pass:");
        userField = new JTextField();
        passwordField = new JPasswordField();
        passPanel.add(message);
        passPanel.add(location);
        passPanel.add(username);
        passPanel.add(userField);
        passPanel.add(password);
        passPanel.add(passwordField);
    }


    public boolean showDialog() {
        int a = JOptionPane.showConfirmDialog(IGV.getMainFrame(), passPanel, "Authentication Required", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);

        if (a == JOptionPane.CANCEL_OPTION) {
            return false;
        } else {
            return true;
        }
    }

    public String getUser() {
        return userField.getText();
    }

    public String getPassword() {
        return new String(passwordField.getPassword());
    }
}
