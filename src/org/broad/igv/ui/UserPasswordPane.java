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
